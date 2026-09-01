"""Decode FlexRay frames from the electrical signal captured with an oscilloscope.

No FlexRay library. The whole protocol is implemented here, from the voltage
threshold to both CRC polynomials.

    python flexray_decode.py capture.csv
    python flexray_decode.py capture.csv --plot

FlexRay does not use bit stuffing the way CAN does. Instead it sends two service
bits, high then low, in front of EVERY byte: the Byte Start Sequence. Its
falling edge is a synchronisation point every 10 bits, and a real receiver uses
it to realign. So does this decoder, which is what makes decoding possible with
only 8 samples per bit.
"""

import argparse

import numpy as np
import pandas as pd

THRESHOLD = 0.30               # V on the difference between the two wires

HEADER_CRC_POLYNOMIAL = 0x385  # x11 + x9 + x8 + x7 + x2 + 1
HEADER_CRC_SEED = 0x01A

FRAME_CRC_POLYNOMIAL = 0x5D6DCB
CHANNEL_SEED = {"A": 0xFEDCBA, "B": 0xABCDEF}


# ------------------------------------------------------------- physical layer

def load(path):
    """Read a PicoScope CSV: time in us, both wires in volts."""
    table = pd.read_csv(path, skiprows=[1]).dropna()
    t = table.iloc[:, 0].to_numpy(dtype=float)
    bp = table.iloc[:, 1].to_numpy(dtype=float)
    bm = table.iloc[:, 2].to_numpy(dtype=float)
    return t, bp - bm


def frame_bounds(difference):
    """Where the frame starts and ends.

    While idle, both wires sit at the same potential and their difference is
    close to zero. The frame is the stretch where that difference departs from
    zero.
    """
    active = np.flatnonzero(np.abs(difference) > THRESHOLD)
    return int(active[0]), int(active[-1])


def to_levels(difference):
    """Voltage difference to logic bits, with hysteresis.

    On this capture a positive difference corresponds to the low level.
    """
    state = 0 if difference[0] > 0 else 1
    out = np.zeros(len(difference), dtype=int)

    for i, v in enumerate(difference):
        if v > THRESHOLD:
            state = 0
        elif v < -THRESHOLD:
            state = 1
        out[i] = state

    return out


def edges(levels, t):
    return t[np.flatnonzero(np.diff(levels)) + 1]


def measure_bit_time(gaps):
    """The bit time, from the gaps between edges.

    With 8 samples per bit each edge is located to within one sample, which is
    12% of a bit. So the shortest gap is no good: it can come out short by pure
    rounding. The most *frequent* gap is, because the error falls on both sides
    and the central value wins.
    """
    values, counts = np.unique(np.round(gaps, 6), return_counts=True)
    bit = float(values[np.argmax(counts)])

    for _ in range(6):                      # every gap lasts a whole number
        multiples = np.round(gaps / bit)    # of bits
        multiples[multiples < 1] = 1
        fits = np.abs(gaps / bit - multiples) < 0.30
        bit = gaps[fits].sum() / multiples[fits].sum()

    return float(bit)


def level_at(t, levels, instant):
    if instant < t[0] or instant > t[-1]:
        return None
    return int(levels[int(np.searchsorted(t, instant))])


def nearest_edge(edge_times, expected, tolerance):
    """The real edge closest to the expected one, if it falls within tolerance."""
    if len(edge_times) == 0:
        return None
    i = int(np.argmin(np.abs(edge_times - expected)))
    return edge_times[i] if abs(edge_times[i] - expected) <= tolerance else None


def extract_bytes(t, levels, edge_times, bit, start):
    """Walk the frame byte by byte, resynchronising on every BSS.

    Ends cleanly on the Frame End Sequence: a low bit followed by a high one
    where a BSS would be due.
    """
    tss_end = edge_times[0]
    tss_bits = round((tss_end - start) / bit)

    border = tss_end + bit                  # after the TSS comes the FSS, then BSS
    out, warnings, closed = [], [], False

    while True:
        first = level_at(t, levels, border + 0.5 * bit)
        second = level_at(t, levels, border + 1.5 * bit)

        if first is None or second is None:
            warnings.append("the capture ends inside the frame")
            break

        if first == 0 and second == 1:      # frame end sequence
            closed = True
            break

        if first != 1 or second != 0:
            warnings.append(f"byte sequence broken at {border:.2f} us")
            break

        falling = nearest_edge(edge_times, border + bit, 0.5 * bit)
        if falling is None:
            warnings.append(f"synchronisation edge lost at {border:.2f} us")
            break

        octet = [level_at(t, levels, falling + (k + 1.5) * bit) for k in range(8)]
        if None in octet:
            warnings.append("the capture ends halfway through a byte")
            break

        out.append(int("".join(str(b) for b in octet), 2))
        border = falling + 9 * bit

    return tss_bits, out, closed, warnings


# ------------------------------------------------------------- protocol

def to_bits(octets):
    return [(b >> i) & 1 for b in octets for i in range(7, -1, -1)]


def to_int(bits):
    return int("".join(str(b) for b in bits), 2)


def crc(bits, polynomial, seed, width):
    """Binary division: the CRC is the remainder."""
    register = seed
    mask = (1 << width) - 1

    for b in bits:
        feedback = b ^ ((register >> (width - 1)) & 1)
        register = (register << 1) & mask
        if feedback:
            register ^= polynomial

    return register


def decode(octets):
    """Interpret the 40 header bits and verify both CRCs."""
    if len(octets) < 8:
        return None

    header = to_bits(octets[:5])

    identifier = to_int(header[5:16])
    words = to_int(header[16:23])
    header_crc_read = to_int(header[23:34])

    payload = octets[5:5 + words * 2]
    tail = octets[5 + words * 2:5 + words * 2 + 3]

    # Header CRC: 11 bits over sync, startup, identifier and length
    header_crc_input = header[3:5] + header[5:16] + header[16:23]
    header_crc_computed = crc(header_crc_input, HEADER_CRC_POLYNOMIAL,
                              HEADER_CRC_SEED, 11)

    # Frame CRC: 24 bits over the whole header plus the payload
    frame_crc_read = to_int(to_bits(tail)) if len(tail) == 3 else None
    frame_crc_input = header + to_bits(payload)

    channel, frame_crc_computed = None, None
    for name, seed in CHANNEL_SEED.items():
        value = crc(frame_crc_input, FRAME_CRC_POLYNOMIAL, seed, 24)
        if value == frame_crc_read:
            channel, frame_crc_computed = name, value
            break
    if frame_crc_computed is None:
        frame_crc_computed = crc(frame_crc_input, FRAME_CRC_POLYNOMIAL,
                                 CHANNEL_SEED["A"], 24)

    return {
        "reserved": header[0],
        "preamble": header[1],
        "not_null": header[2],
        "sync": header[3],
        "startup": header[4],
        "id": identifier,
        "words": words,
        "payload_bytes": words * 2,
        "payload": payload,
        "cycle": to_int(header[34:40]),
        "header_crc_read": header_crc_read,
        "header_crc_computed": header_crc_computed,
        "header_crc_ok": header_crc_read == header_crc_computed,
        "frame_crc_read": frame_crc_read,
        "frame_crc_computed": frame_crc_computed,
        "frame_crc_ok": frame_crc_read == frame_crc_computed,
        "channel": channel,
    }


def analyse(path):
    t, difference = load(path)
    levels = to_levels(difference)

    a, b = frame_bounds(difference)
    edge_times = edges(levels[a:b + 1], t[a:b + 1])
    bit = measure_bit_time(np.diff(edge_times))

    tss_bits, octets, closed, warnings = extract_bytes(
        t, levels, edge_times, bit, t[a])

    return {
        "samples": len(t),
        "step_us": float(t[1] - t[0]),
        "bit_us": bit,
        "mbits": 1.0 / bit,
        "samples_per_bit": bit / (t[1] - t[0]),
        "tss_bits": tss_bits,
        "bytes_read": octets,
        "closed": closed,
        "warnings": warnings,
        "frame": decode(octets),
    }


# ------------------------------------------------------------- output

def report(r, path):
    print(f"\n{path}")
    print(f"  samples             {r['samples']}")
    print(f"  sampling step       {r['step_us']:.4f} us   ({1 / r['step_us']:.0f} MS/s)")
    print(f"  bit time            {r['bit_us'] * 1000:.1f} ns")
    print(f"  bus speed           {r['mbits']:.2f} Mbit/s")
    print(f"  samples per bit     {r['samples_per_bit']:.1f}")

    print(f"\n  TSS                 {r['tss_bits']} bits")
    print(f"  bytes read          {len(r['bytes_read'])}")
    print(f"  frame end           {'FES found' if r['closed'] else 'NOT FOUND'}")

    for warning in r["warnings"]:
        print(f"  warning: {warning}")

    m = r["frame"]
    if m is None:
        print("\n  incomplete frame: cannot be interpreted")
        return

    kind = []
    if m["sync"]:
        kind.append("sync")
    if m["startup"]:
        kind.append("startup")
    if not m["not_null"]:
        kind.append("null")

    print(f"\n  identifier          {m['id']}")
    print(f"  type                {', '.join(kind) if kind else 'data'}")
    print(f"  cycle counter       {m['cycle']}")
    print(f"  payload             {m['words']} words = {m['payload_bytes']} bytes")

    for i in range(0, len(m["payload"]), 16):
        print(f"     {' '.join(f'{b:02X}' for b in m['payload'][i:i + 16])}")

    print(f"\n  header CRC          0x{m['header_crc_read']:03X} read / "
          f"0x{m['header_crc_computed']:03X} computed   "
          f"{'MATCH' if m['header_crc_ok'] else 'MISMATCH'}")
    print(f"  frame CRC           0x{m['frame_crc_read']:06X} read / "
          f"0x{m['frame_crc_computed']:06X} computed   "
          f"{'MATCH' if m['frame_crc_ok'] else 'MISMATCH'}")

    if m["channel"]:
        print(f"  channel             {m['channel']}   "
              f"(inferred from the CRC seed that matches)")


def plot(path, out="flexray.png"):
    import matplotlib.pyplot as plt

    table = pd.read_csv(path, skiprows=[1]).dropna()
    t = table.iloc[:, 0].to_numpy(dtype=float)
    bp = table.iloc[:, 1].to_numpy(dtype=float)
    bm = table.iloc[:, 2].to_numpy(dtype=float)

    fig, (top, bottom) = plt.subplots(2, 1, figsize=(13, 6), sharex=True)

    top.plot(t, bp, linewidth=0.7, label="BP")
    top.plot(t, bm, linewidth=0.7, label="BM")
    top.set_ylabel("V")
    top.legend(loc="upper right")
    top.set_title("Both wires")

    bottom.plot(t, bp - bm, linewidth=0.7, color="black")
    bottom.axhline(THRESHOLD, linestyle="--", linewidth=0.8)
    bottom.axhline(-THRESHOLD, linestyle="--", linewidth=0.8)
    bottom.set_ylabel("V")
    bottom.set_xlabel("time (us)")
    bottom.set_title("BP minus BM, with both thresholds")

    for axis in (top, bottom):
        axis.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(out, dpi=110)
    print(f"\n  plot saved to {out}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", help="CSV exported from PicoScope")
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()

    report(analyse(args.csv), args.csv)

    if args.plot:
        plot(args.csv)


if __name__ == "__main__":
    main()
