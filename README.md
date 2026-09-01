# flexray-physical-decoder

Decodes FlexRay frames from the electrical signal captured with an oscilloscope,
and generates synthetic FlexRay signals to test it. No FlexRay library: the whole
protocol is implemented from scratch, from the voltage threshold to both CRC
polynomials.

![tests](https://github.com/andriusvillena-cmd/flexray-physical-decoder/actions/workflows/tests.yml/badge.svg)

---

## What it does

Input: a PicoScope CSV with three columns — time and the two wires of the
differential pair, in volts.

```
python flexray_decode.py my_capture.csv
```

```
  sampling step       0.0125 us   (80 MS/s)
  bit time            100.0 ns
  bus speed           10.00 Mbit/s
  samples per bit     8.0

  TSS                 11 bits
  bytes read          40
  frame end           FES found

  identifier          1
  type                sync, startup
  cycle counter       45
  payload             16 words = 32 bytes

  header CRC          0x0F2 read / 0x0F2 computed   MATCH
  frame CRC           0xA011D3 read / 0xA011D3 computed   MATCH
  channel             A   (inferred from the CRC seed that matches)
```

And the other way round:

```
python flexray_gen.py --id 42 --cycle 7 --bytes 16 --pattern count
```

---

## Why FlexRay is not decoded like CAN

**There is no bit stuffing.** CAN forbids more than 5 identical bits in a row
and inserts an opposite one to force an edge. FlexRay does something else: in
front of **every byte** it sends two service bits, one high and one low, called
the Byte Start Sequence. The falling edge of those two bits is a synchronisation
point every 10 bits.

That changes the design of the decoder. In CAN you measure the bit time once and
sample the whole frame counting from the start bit. Do that here and drift eats
you: at 10 Mbit/s with 8 samples per bit, a 1% error in the bit time shifts you
by more than half a bit before the end of the frame. So this decoder
**realigns on every BSS**, exactly as a real receiver does, and that is why it
still works with 4 samples per bit.

**There are two CRCs, not one.** An 11-bit one covering only the header, and a
24-bit one covering header and payload. The header is protected separately
because a node has to be able to trust the identifier and the length *before* it
has received the whole frame.

**And the frame CRC does not start from zero.** It starts from a different
constant depending on the physical channel: one for A and another for B. That
has a neat practical consequence: try both, and if only one matches you have
worked out which channel the capture came from without anyone telling you.

---

## How it works, step by step

**1. Differential.** The two wires are subtracted. While idle they sit at the
same potential and the difference is close to zero, which makes it possible to
find where the frame starts and ends without looking for patterns.

**2. Threshold with hysteresis.** Two lines at ±0.3 V, and between them the
previous state is held.

**3. Bit time.** Measured from the signal itself, and **not** from the shortest
gap between edges. With 8 samples per bit each edge is located to within one
sample, which is 12% of a bit, and the minimum comes out short by pure rounding.
The most *frequent* gap is used instead, where the error falls on both sides,
and it is then refined by requiring every gap to last a whole number of bits.

**4. TSS.** The frame opens with a stretch of low level, between 3 and 15 bits.

**5. Byte-by-byte walk.** FSS, and then for each byte: the BSS is checked, its
real falling edge is located, the decoder anchors there and reads the 8 bits at
their centres. That anchoring is what prevents drift.

**6. FES.** The frame closes with a low bit followed by a high one where a BSS
would be due. The decoder tells that apart from a synchronisation error.

**7. Header.** 40 bits: reserved, payload preamble indicator, null frame
indicator, sync, startup, 11-bit identifier, 7-bit payload length in two-byte
words, 11-bit header CRC, and a 6-bit cycle counter.

**8. Both CRCs.** Recomputed and compared with the transmitted ones. The header
one with polynomial 0x385 from 0x01A, over 20 bits. The frame one with
polynomial 0x5D6DCB, over the whole header plus the payload, trying both channel
seeds.

---

## The data

**This repository contains no oscilloscope capture.**

The decoder was developed and validated against a real capture from the
PicoScope waveform library: a startup frame, channel A, identifier 1, cycle 45,
32 payload bytes all zero, taken at 80 MS/s. The reuse terms for that library
are not published, so the file is not redistributed.

What is published is `flexray_gen.py`, which produces equivalent signals. It
allows the decoder to be tested on any machine, and under conditions a single
capture cannot provide: different identifiers, payload lengths from 2 to 254
bytes, both channels, lower sample rates, noise and reduced amplitude.

Writing the generator was not a detour. Building a frame forces you to know what
goes in every bit, and leaves the protocol implemented from both sides.

---

## The closed circle, and how it is broken

Testing a decoder with data from your own generator has an obvious risk: if the
two share the same misunderstanding, the tests go green and the result is wrong.

That is why there are two **anchor** tests comparing against values measured on
the real capture:

- The 5 header bytes the generator produces are `38 01 20 3C AD`.
- Both CRCs come out as `0x0F2` and `0xA011D3`.

Those numbers came off the oscilloscope before the generator existed. If either
module drifts from the standard, those tests go red even while everything else
keeps agreeing with itself.

---

## Tests

```
pytest -v
```

Eighteen tests. Green means the code works.

| Group | What it checks |
|---|---|
| Anchor | Header and both CRCs equal to the ones measured on the real capture |
| Physical layer | 10 Mbit/s, TSS of 11 and of 5 bits, FES, hysteresis, bit time measurement against an anomalous gap |
| Robustness | 4 samples per bit, increased noise, reduced amplitude |
| Round trip | Identifiers 1, 42 and 2047; cycles 0, 7 and 63; payloads of 2, 16 and 254 bytes |
| Channel | The channel is inferred from the CRC seed that matches |
| Fault detection | A corrupted bit in the payload breaks the frame CRC but not the header CRC; one in the header breaks its own |

The four-samples-per-bit test is the one that says the most: at 40 MS/s each bit
gets four points, and the only reason the frame still comes out whole is the
resynchronisation on every BSS.

---

## Installation

```
pip install -r requirements.txt
```

Python 3.10 or later. numpy, pandas, matplotlib and pytest.

---

## Related

[can-physical-decoder](https://github.com/andriusvillena-cmd/can-physical-decoder) —
the same for CAN: from the oscilloscope signal to the message, with CRC
verification and detection of frames nobody acknowledged.

## Reference standards

ISO 17458-2 (data link layer) · ISO 17458-4 (physical layer)
