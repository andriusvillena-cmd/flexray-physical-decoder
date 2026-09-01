"""Genera una senal de FlexRay sintetica, en el mismo formato que PicoScope.

Sirve para dos cosas. Una: probar el decodificador sin depender de ninguna
captura ajena, que es lo que permite que las pruebas corran solas en cualquier
maquina. Y dos: entender el protocolo desde el otro lado, porque construir una
trama obliga a saber exactamente que va en cada bit.

    python flexray_gen.py                       una trama de arranque
    python flexray_gen.py --salida mia.csv --id 42 --ciclo 7 --bytes 16
"""

import argparse

import numpy as np

TENSION_REPOSO = 2.5           # V, los dos hilos juntos cuando no se transmite

CRC_CABECERA = 0x385
INICIO_CABECERA = 0x01A

CRC_TRAMA = 0x5D6DCB
INICIO_CANAL = {"A": 0xFEDCBA, "B": 0xABCDEF}


def crc(bits, polinomio, inicio, ancho):
    registro = inicio
    mascara = (1 << ancho) - 1

    for b in bits:
        siguiente = b ^ ((registro >> (ancho - 1)) & 1)
        registro = (registro << 1) & mascara
        if siguiente:
            registro ^= polinomio

    return registro


def a_bits(valor, ancho):
    return [(valor >> i) & 1 for i in range(ancho - 1, -1, -1)]


def bits_de_bytes(octetos):
    return [(b >> i) & 1 for b in octetos for i in range(7, -1, -1)]


def bytes_de_bits(bits):
    return [int("".join(str(b) for b in bits[i:i + 8]), 2)
            for i in range(0, len(bits), 8)]


# ------------------------------------------------------------- la trama

def construir_trama(ident, ciclo, carga, sincronismo=1, arranque=1, canal="A"):
    """Devuelve los 40 bits de cabecera, la carga util y los 3 bytes de CRC."""
    if len(carga) % 2:
        raise ValueError("la carga util va en palabras de 2 bytes")

    palabras = len(carga) // 2

    # Los 20 bits que protege el CRC de cabecera
    entrada = ([sincronismo, arranque]
               + a_bits(ident, 11)
               + a_bits(palabras, 7))
    crc_cab = crc(entrada, CRC_CABECERA, INICIO_CABECERA, 11)

    cabecera = ([0]                      # reservado
                + [0]                    # indicador de preambulo
                + [1]                    # no es trama nula
                + [sincronismo]
                + [arranque]
                + a_bits(ident, 11)
                + a_bits(palabras, 7)
                + a_bits(crc_cab, 11)
                + a_bits(ciclo, 6))

    crc_tr = crc(cabecera + bits_de_bytes(carga),
                 CRC_TRAMA, INICIO_CANAL[canal], 24)

    return cabecera, list(carga), a_bits(crc_tr, 24)


def cadena_de_bits(cabecera, carga, cola_bits, bits_tss=11):
    """Monta la secuencia completa que viaja por el cable.

    TSS, FSS, y luego cada byte precedido por su BSS (alto y bajo).
    Cierra con el FES: bajo y alto.
    """
    octetos = bytes_de_bits(cabecera) + carga + bytes_de_bits(cola_bits)

    cadena = [0] * bits_tss              # TSS: nivel bajo mantenido
    cadena += [1]                        # FSS

    for b in octetos:
        cadena += [1, 0]                 # BSS
        cadena += [(b >> i) & 1 for i in range(7, -1, -1)]

    cadena += [0, 1]                     # FES
    return cadena


# ------------------------------------------------------------- la senal

def a_tension(cadena, bit_us=0.1, muestreo_ms=80.0, amplitud=0.9,
              reposo_us=5.0, subida_ns=25.0, ruido=0.006, semilla=0):
    """Convierte la cadena de bits en dos tensiones, como las veria una sonda.

    Nivel bajo -> el primer hilo por encima del segundo.
    Los flancos no son verticales: se suavizan para imitar el tiempo de subida.
    """
    paso = 1.0 / muestreo_ms                     # us entre muestras
    por_bit = int(round(bit_us / paso))
    if por_bit < 3:
        raise ValueError("muestreo demasiado bajo para este tiempo de bit")

    huecos = int(round(reposo_us / paso))

    nivel = np.concatenate([
        np.zeros(huecos),
        np.repeat([1.0 if b == 0 else -1.0 for b in cadena], por_bit),
        np.zeros(huecos),
    ])

    # Suavizado de flancos: media movil de la anchura del tiempo de subida
    ancho = max(1, int(round((subida_ns / 1000.0) / paso)))
    if ancho > 1:
        nucleo = np.ones(ancho) / ancho
        nivel = np.convolve(nivel, nucleo, mode="same")

    generador = np.random.default_rng(semilla)
    t = (np.arange(len(nivel)) * paso) - reposo_us

    canal_a = TENSION_REPOSO + nivel * amplitud / 2
    canal_b = TENSION_REPOSO - nivel * amplitud / 2

    canal_a += generador.normal(0, ruido, len(nivel))
    canal_b += generador.normal(0, ruido, len(nivel))

    return t, canal_a, canal_b


def generar(ident=1, ciclo=45, carga=None, canal="A", sincronismo=1,
            arranque=1, bits_tss=11, **opciones):
    """Atajo: de los parametros de la trama a las tres columnas de la senal."""
    if carga is None:
        carga = [0x00] * 32

    cabecera, carga, cola = construir_trama(
        ident, ciclo, carga, sincronismo, arranque, canal)

    return a_tension(cadena_de_bits(cabecera, carga, cola, bits_tss), **opciones)


def escribir_csv(ruta, t, canal_a, canal_b):
    """Mismo formato que exporta PicoScope 6: dos lineas de cabecera y una vacia."""
    with open(ruta, "w", encoding="utf-8") as f:
        f.write("Tiempo,Canal A,Canal C\n(us),(V),(V)\n\n")
        for i in range(len(t)):
            f.write(f"{t[i]:.8f},{canal_a[i]:.8f},{canal_b[i]:.8f}\n")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--salida", default="flexray_sintetica.csv")
    p.add_argument("--id", type=int, default=1, dest="ident")
    p.add_argument("--ciclo", type=int, default=45)
    p.add_argument("--bytes", type=int, default=32, dest="n_bytes")
    p.add_argument("--canal", choices=["A", "B"], default="A")
    p.add_argument("--muestreo", type=float, default=80.0, help="MS/s")
    p.add_argument("--patron", choices=["ceros", "cuenta", "azar"], default="ceros")
    args = p.parse_args()

    if args.patron == "ceros":
        carga = [0x00] * args.n_bytes
    elif args.patron == "cuenta":
        carga = [i & 0xFF for i in range(args.n_bytes)]
    else:
        carga = list(np.random.default_rng(0).integers(0, 256, args.n_bytes))

    t, a, b = generar(ident=args.ident, ciclo=args.ciclo, carga=carga,
                      canal=args.canal, muestreo_ms=args.muestreo)

    escribir_csv(args.salida, t, a, b)

    print(f"{args.salida}")
    print(f"  ID {args.ident}, ciclo {args.ciclo}, {len(carga)} bytes, canal {args.canal}")
    print(f"  {len(t)} muestras a {args.muestreo:g} MS/s")


if __name__ == "__main__":
    main()
