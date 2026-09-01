"""Regression tests for the FlexRay decoder.

Every signal is produced by the generator, so the tests run on any machine with
no external capture file.

To stop that being a closed circle -- generator and decoder agreeing on the same
misunderstanding -- two anchor tests compare against values measured on a real
oscilloscope capture, taken separately and not included in this repository.

    pytest -v
"""

import numpy as np
import pytest

import flexray_decode as fd
import flexray_gen as fg


# Measured on a real PicoScope capture, channel A: startup frame, identifier 1,
# cycle 45, 32 payload bytes all zero.
REAL_HEADER = [0x38, 0x01, 0x20, 0x3C, 0xAD]
REAL_HEADER_CRC = 0x0F2
REAL_FRAME_CRC = 0xA011D3


def signal(tmp_path, name="signal.csv", **options):
    """Generate a signal and leave it in a temporary CSV."""
    path = tmp_path / name
    t, a, b = fg.generate(**options)
    fg.write_csv(path, t, a, b)
    return str(path)


# ------------------------------------------------------- anchored in reality

def test_header_matches_the_real_capture():
    """The 5 header bytes the generator produces are the ones measured."""
    header, _, _ = fg.build_frame(1, 45, [0x00] * 32)
    assert fg.bytes_from_bits(header) == REAL_HEADER


def test_both_crcs_match_the_real_capture(tmp_path):
    m = fd.analyse(signal(tmp_path))["frame"]

    assert m["header_crc_read"] == REAL_HEADER_CRC
    assert m["frame_crc_read"] == REAL_FRAME_CRC
    assert m["header_crc_ok"]
    assert m["frame_crc_ok"]


# ------------------------------------------------------- physical layer

def test_speed_and_sampling(tmp_path):
    r = fd.analyse(signal(tmp_path))

    assert r["bit_us"] == pytest.approx(0.1, abs=0.001)
    assert r["mbits"] == pytest.approx(10.0, abs=0.1)
    assert r["samples_per_bit"] == pytest.approx(8.0, abs=0.1)


def test_tss_and_frame_end(tmp_path):
    r = fd.analyse(signal(tmp_path, tss_bits=11))

    assert r["tss_bits"] == 11
    assert r["closed"] is True
    assert r["warnings"] == []


def test_a_different_tss_length(tmp_path):
    """The standard allows between 3 and 15 start bits."""
    assert fd.analyse(signal(tmp_path, tss_bits=5))["tss_bits"] == 5


def test_hysteresis_holds_the_state_in_the_dead_band():
    trace = np.array([1.0, -1.0, 0.0, 0.0, 1.0, 0.0, -1.0])
    assert list(fd.to_levels(trace)) == [0, 1, 1, 1, 0, 0, 1]


def test_bit_time_ignores_an_anomalous_gap():
    gaps = np.array([0.1] * 30 + [0.2] * 8 + [0.037])
    assert fd.measure_bit_time(gaps) == pytest.approx(0.1, abs=0.002)


def test_survives_only_four_samples_per_bit(tmp_path):
    """At 40 MS/s each bit gets four samples. Resynchronising on every BSS is
    the only reason the frame still comes out whole."""
    r = fd.analyse(signal(tmp_path, sample_rate_ms=40.0))

    assert r["samples_per_bit"] == pytest.approx(4.0, abs=0.1)
    assert r["frame"]["frame_crc_ok"]


def test_survives_noise(tmp_path):
    r = fd.analyse(signal(tmp_path, noise=0.05, seed=7))
    assert r["frame"]["frame_crc_ok"]


def test_survives_low_amplitude(tmp_path):
    r = fd.analyse(signal(tmp_path, amplitude=0.7))
    assert r["frame"]["frame_crc_ok"]


# ------------------------------------------------------- round trip

@pytest.mark.parametrize("identifier,cycle,n", [(1, 0, 2), (42, 7, 16), (2047, 63, 254)])
def test_round_trip(tmp_path, identifier, cycle, n):
    payload = [(i * 7 + 3) & 0xFF for i in range(n)]
    path = signal(tmp_path, identifier=identifier, cycle=cycle, payload=payload)

    m = fd.analyse(path)["frame"]

    assert m["id"] == identifier
    assert m["cycle"] == cycle
    assert m["payload_bytes"] == n
    assert m["payload"] == payload
    assert m["header_crc_ok"]
    assert m["frame_crc_ok"]


def test_frame_indicators(tmp_path):
    m = fd.analyse(signal(tmp_path, sync=0, startup=0))["frame"]

    assert m["sync"] == 0
    assert m["startup"] == 0
    assert m["reserved"] == 0
    assert m["not_null"] == 1


def test_channel_is_inferred_from_the_crc_seed(tmp_path):
    """The frame CRC starts from a different constant on each channel."""
    assert fd.analyse(signal(tmp_path, channel="A"))["frame"]["channel"] == "A"
    assert fd.analyse(signal(tmp_path, channel="B", name="b.csv"))["frame"]["channel"] == "B"


# ------------------------------------------------------- fault detection

def test_a_corrupted_byte_breaks_the_frame_crc():
    header, payload, tail = fg.build_frame(1, 45, [0x11] * 8)
    payload[3] ^= 0x01                                 # one bit of the payload

    octets = fg.bytes_from_bits(header) + payload + fg.bytes_from_bits(tail)

    m = fd.decode(octets)
    assert m["header_crc_ok"]                          # the header is intact
    assert not m["frame_crc_ok"]                       # the data is not


def test_a_corrupted_header_breaks_its_own_crc():
    header, payload, tail = fg.build_frame(1, 45, [0x11] * 8)
    octets = fg.bytes_from_bits(header) + payload + fg.bytes_from_bits(tail)
    octets[1] ^= 0x01                                  # one bit of the identifier

    assert not fd.decode(octets)["header_crc_ok"]


def test_too_short_a_frame_is_not_interpreted():
    assert fd.decode([0x38, 0x01, 0x20]) is None
