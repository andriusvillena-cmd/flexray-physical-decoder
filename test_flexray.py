"""Pruebas de regresion del decodificador de FlexRay.

Todas las senales se fabrican con el generador, asi que las pruebas corren en
cualquier maquina sin depender de ninguna captura externa.

Para que eso no se convierta en un circulo cerrado -- generador y decodificador
poniendose de acuerdo en el mismo error -- hay dos pruebas de anclaje que
comparan contra valores leidos de una captura real de osciloscopio, medida
aparte y no incluida en el repositorio.

    pytest -v
"""

import numpy as np
import pytest

import flexray_decode as fd
import flexray_gen as fg


# Valores medidos en una captura real de PicoScope, canal A:
# trama de arranque, identificador 1, ciclo 45, 32 bytes de carga a cero.
CABECERA_REAL = [0x38, 0x01, 0x20, 0x3C, 0xAD]
CRC_CABECERA_REAL = 0x0F2
CRC_TRAMA_REAL = 0xA011D3


def archivo(tmp_path, nombre="senal.csv", **opciones):
    """Genera una senal y la deja en un CSV temporal."""
    ruta = tmp_path / nombre
    t, a, b = fg.generar(**opciones)
    fg.escribir_csv(ruta, t, a, b)
    return str(ruta)


# ------------------------------------------------------- anclaje en lo real

def test_cabecera_igual_a_la_captura_real():
    """Los 5 bytes de cabecera que fabrica el generador son los medidos."""
    cabecera, _, _ = fg.construir_trama(1, 45, [0x00] * 32)
    assert fg.bytes_de_bits(cabecera) == CABECERA_REAL


def test_los_dos_crc_igual_que_en_la_captura_real(tmp_path):
    m = fd.analizar(archivo(tmp_path))["trama"]

    assert m["crc_cab_leido"] == CRC_CABECERA_REAL
    assert m["crc_trama_leido"] == CRC_TRAMA_REAL
    assert m["crc_cab_ok"]
    assert m["crc_trama_ok"]


# ------------------------------------------------------- capa fisica

def test_velocidad_y_muestreo(tmp_path):
    r = fd.analizar(archivo(tmp_path))

    assert r["bit_us"] == pytest.approx(0.1, abs=0.001)
    assert r["mbits"] == pytest.approx(10.0, abs=0.1)
    assert r["muestras_por_bit"] == pytest.approx(8.0, abs=0.1)


def test_tss_y_fin_de_trama(tmp_path):
    r = fd.analizar(archivo(tmp_path, bits_tss=11))

    assert r["bits_tss"] == 11
    assert r["cerrada"] is True
    assert r["avisos"] == []


def test_tss_de_otra_longitud(tmp_path):
    """La norma admite entre 3 y 15 bits de arranque."""
    assert fd.analizar(archivo(tmp_path, bits_tss=5))["bits_tss"] == 5


def test_histeresis_conserva_el_estado_en_la_zona_muerta():
    senal = np.array([1.0, -1.0, 0.0, 0.0, 1.0, 0.0, -1.0])
    assert list(fd.digitalizar(senal)) == [0, 1, 1, 1, 0, 0, 1]


def test_medir_bit_ignora_un_hueco_anomalo():
    huecos = np.array([0.1] * 30 + [0.2] * 8 + [0.037])
    assert fd.medir_bit(huecos) == pytest.approx(0.1, abs=0.002)


def test_aguanta_con_solo_cuatro_muestras_por_bit(tmp_path):
    """A 40 MS/s cada bit tiene 4 muestras. La sincronizacion en cada BSS
    es justo lo que permite que siga saliendo."""
    r = fd.analizar(archivo(tmp_path, muestreo_ms=40.0))

    assert r["muestras_por_bit"] == pytest.approx(4.0, abs=0.1)
    assert r["trama"]["crc_trama_ok"]


def test_aguanta_ruido(tmp_path):
    r = fd.analizar(archivo(tmp_path, ruido=0.05, semilla=7))
    assert r["trama"]["crc_trama_ok"]


def test_aguanta_amplitud_baja(tmp_path):
    r = fd.analizar(archivo(tmp_path, amplitud=0.7))
    assert r["trama"]["crc_trama_ok"]


# ------------------------------------------------------- ida y vuelta

@pytest.mark.parametrize("ident,ciclo,n", [(1, 0, 2), (42, 7, 16), (2047, 63, 254)])
def test_ida_y_vuelta(tmp_path, ident, ciclo, n):
    carga = [(i * 7 + 3) & 0xFF for i in range(n)]
    ruta = archivo(tmp_path, ident=ident, ciclo=ciclo, carga=carga)

    m = fd.analizar(ruta)["trama"]

    assert m["id"] == ident
    assert m["ciclo"] == ciclo
    assert m["bytes_carga"] == n
    assert m["carga"] == carga
    assert m["crc_cab_ok"]
    assert m["crc_trama_ok"]


def test_indicadores_de_trama(tmp_path):
    m = fd.analizar(archivo(tmp_path, sincronismo=0, arranque=0))["trama"]

    assert m["sincronismo"] == 0
    assert m["arranque"] == 0
    assert m["reservado"] == 0
    assert m["no_nula"] == 1


def test_el_canal_se_deduce_de_la_semilla_del_crc(tmp_path):
    """El CRC de trama arranca de una constante distinta en cada canal."""
    assert fd.analizar(archivo(tmp_path, canal="A"))["trama"]["canal"] == "A"
    assert fd.analizar(archivo(tmp_path, canal="B", nombre="b.csv"))["trama"]["canal"] == "B"


# ------------------------------------------------------- deteccion de fallos

def test_un_byte_corrompido_rompe_el_crc_de_trama():
    cabecera, carga, cola = fg.construir_trama(1, 45, [0x11] * 8)
    carga[3] ^= 0x01                                   # un bit de la carga util

    cadena = fg.cadena_de_bits(cabecera, carga, cola)
    octetos = (fg.bytes_de_bits(cabecera) + carga + fg.bytes_de_bits(cola))

    m = fd.decodificar(octetos)
    assert m["crc_cab_ok"]                             # la cabecera esta intacta
    assert not m["crc_trama_ok"]                       # los datos no


def test_una_cabecera_corrompida_rompe_su_propio_crc():
    cabecera, carga, cola = fg.construir_trama(1, 45, [0x11] * 8)
    octetos = fg.bytes_de_bits(cabecera) + carga + fg.bytes_de_bits(cola)
    octetos[1] ^= 0x01                                 # un bit del identificador

    assert not fd.decodificar(octetos)["crc_cab_ok"]


def test_trama_demasiado_corta_no_se_interpreta():
    assert fd.decodificar([0x38, 0x01, 0x20]) is None
