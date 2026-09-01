"""Generate a synthetic FlexRay signal, in the same format PicoScope exports.

It serves two purposes. One: testing the decoder without depending on anyone
else's capture, which is what lets the test suite run on any machine. Two:
understanding the protocol from the other side, because building a frame forces
you to know exactly what goes in every bit.

    python flexray_gen.py
    python flexray_gen.py --out mine.csv --id 42 --cycle 7 --bytes 16
"""

import argparse

import numpy as np

IDLE_VOLTAGE = 2.5             # V, both wires together while nothing transmits

HEADER_CRC_POLYNOMIAL = 0x385
HEADER_CRC_SEED = 0x01A

FRAME_CRC_POLYNOMIAL = 0x5D6DCB
CHANNEL_SEED = {"A": 0xFEDCBA, "B": 0xABCDEF}


def crc(bits, polynomial, seed, width):
    register = seed
    mask = (1 << width) - 1

    for b in bits:
        feedback = b ^ ((register >> (width - 1)) & 1)
        register = (register << 1) & mask
        if feedback:
            register ^= polynomial

    return register


def to_bits(value, width):
    return [(value >> i) & 1 for i in range(width - 1, -1, -1)]


def bits_from_bytes(octets):
    return [(b >> i) & 1 for b in octets for i in range(7, -1, -1)]


def bytes_from_bits(bits):
    return [int("".join(str(b) for b in bits[i:i + 8]), 2)
            for i in range(0, len(bits), 8)]


# ------------------------------------------------------------- the frame

def build_frame(identifier, cycle, payload, sync=1, startup=1, channel="A"):
    """Return the 40 header bits, the payload, and the 3 CRC bytes."""
    if len(payload) % 2:
        raise ValueError("the payload goes in two-byte words")

    words = len(payload) // 2

    # The 20 bits the header CRC protects
    crc_input = ([sync, startup]
                 + to_bits(identifier, 11)
                 + to_bits(words, 7))
    header_crc = crc(crc_input, HEADER_CRC_POLYNOMIAL, HEADER_CRC_SEED, 11)

    header = ([0]                        # reserved
              + [0]                      # payload preamble indicator
              + [1]                      # not a null frame
              + [sync]
              + [startup]
              + to_bits(identifier, 11)
              + to_bits(words, 7)
              + to_bits(header_crc, 11)
              + to_bits(cycle, 6))

    frame_crc = crc(header + bits_from_bytes(payload),
                    FRAME_CRC_POLYNOMIAL, CHANNEL_SEED[channel], 24)

    return header, list(payload), to_bits(frame_crc, 24)


def bit_stream(header, payload, tail_bits, tss_bits=11):
    """Assemble the full sequence that travels down the wire.

    TSS, FSS, then every byte preceded by its BSS (high then low). Closes with
    the FES: low then high.
    """
    octets = bytes_from_bits(header) + payload + bytes_from_bits(tail_bits)

    stream = [0] * tss_bits              # TSS: sustained low level
    stream += [1]                        # FSS

    for b in octets:
        stream += [1, 0]                 # BSS
        stream += [(b >> i) & 1 for i in range(7, -1, -1)]

    stream += [0, 1]                     # FES
    return stream


# ------------------------------------------------------------- the signal

def to_voltage(stream, bit_us=0.1, sample_rate_ms=80.0, amplitude=0.9,
               idle_us=5.0, rise_ns=25.0, noise=0.006, seed=0):
    """Turn the bit stream into two voltages, as a probe would see them.

    Low level -> the first wire above the second one.
    Edges are not vertical: they are smoothed to imitate a real rise time.
    """
    step = 1.0 / sample_rate_ms                  # us between samples
    per_bit = int(round(bit_us / step))
    if per_bit < 3:
        raise ValueError("sample rate too low for this bit time")

    pad = int(round(idle_us / step))

    level = np.concatenate([
        np.zeros(pad),
        np.repeat([1.0 if b == 0 else -1.0 for b in stream], per_bit),
        np.zeros(pad),
    ])

    # Edge smoothing: a moving average as wide as the rise time
    width = max(1, int(round((rise_ns / 1000.0) / step)))
    if width > 1:
        level = np.convolve(level, np.ones(width) / width, mode="same")

    rng = np.random.default_rng(seed)
    t = (np.arange(len(level)) * step) - idle_us

    wire_a = IDLE_VOLTAGE + level * amplitude / 2
    wire_b = IDLE_VOLTAGE - level * amplitude / 2

    wire_a += rng.normal(0, noise, len(level))
    wire_b += rng.normal(0, noise, len(level))

    return t, wire_a, wire_b


def generate(identifier=1, cycle=45, payload=None, channel="A", sync=1,
             startup=1, tss_bits=11, **options):
    """Shortcut: from frame parameters to the three columns of the signal."""
    if payload is None:
        payload = [0x00] * 32

    header, payload, tail = build_frame(
        identifier, cycle, payload, sync, startup, channel)

    return to_voltage(bit_stream(header, payload, tail, tss_bits), **options)


def write_csv(path, t, wire_a, wire_b):
    """The same format PicoScope 6 exports: two header lines and a blank one."""
    with open(path, "w", encoding="utf-8") as f:
        f.write("Time,Channel A,Channel C\n(us),(V),(V)\n\n")
        for i in range(len(t)):
            f.write(f"{t[i]:.8f},{wire_a[i]:.8f},{wire_b[i]:.8f}\n")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", default="flexray_synthetic.csv")
    p.add_argument("--id", type=int, default=1, dest="identifier")
    p.add_argument("--cycle", type=int, default=45)
    p.add_argument("--bytes", type=int, default=32, dest="n_bytes")
    p.add_argument("--channel", choices=["A", "B"], default="A")
    p.add_argument("--sample-rate", type=float, default=80.0, help="MS/s")
    p.add_argument("--pattern", choices=["zeros", "count", "random"], default="zeros")
    args = p.parse_args()

    if args.pattern == "zeros":
        payload = [0x00] * args.n_bytes
    elif args.pattern == "count":
        payload = [i & 0xFF for i in range(args.n_bytes)]
    else:
        payload = list(np.random.default_rng(0).integers(0, 256, args.n_bytes))

    t, a, b = generate(identifier=args.identifier, cycle=args.cycle,
                       payload=payload, channel=args.channel,
                       sample_rate_ms=args.sample_rate)

    write_csv(args.out, t, a, b)

    print(f"{args.out}")
    print(f"  id {args.identifier}, cycle {args.cycle}, {len(payload)} bytes, "
          f"channel {args.channel}")
    print(f"  {len(t)} samples at {args.sample_rate:g} MS/s")


if __name__ == "__main__":
    main()
