"""FlexRay: de la senal del osciloscopio a la trama verificada.

Sin librerias de FlexRay. Todo el protocolo esta implementado aqui, del umbral
de tension a los dos polinomios de CRC.

    python flexray_decode.py FlexRay_Trace.csv
    python flexray_decode.py FlexRay_Trace.csv --grafica

FlexRay no usa relleno de bits como CAN. Delante de CADA byte manda dos bits de
servicio, alto y bajo, llamados BSS. Su flanco de bajada es un punto de
sincronizacion cada 10 bits, y el receptor lo aprovecha para recolocarse. Aqui
se hace lo mismo, que es lo que permite decodificar con solo 8 muestras por bit.
"""

import argparse

import numpy as np
import pandas as pd

UMBRAL = 0.30                  # V en la diferencia entre los dos hilos

CRC_CABECERA = 0x385           # x11+x9+x8+x7+x2+1
INICIO_CABECERA = 0x01A

CRC_TRAMA = 0x5D6DCB           # polinomio de 24 bits
INICIO_CANAL = {"A": 0xFEDCBA, "B": 0xABCDEF}


# ------------------------------------------------------------- capa fisica

def cargar(ruta):
    """Lee el CSV de PicoScope: tiempo en us y los dos hilos en voltios."""
    tabla = pd.read_csv(ruta, skiprows=[1]).dropna()
    t = tabla.iloc[:, 0].to_numpy(dtype=float)
    bp = tabla.iloc[:, 1].to_numpy(dtype=float)
    bm = tabla.iloc[:, 2].to_numpy(dtype=float)
    return t, bp - bm


def acotar(dif):
    """Donde empieza y acaba la trama.

    En reposo los dos hilos estan al mismo potencial y su diferencia es casi
    cero. La trama es el tramo donde esa diferencia se separa de cero.
    """
    activo = np.flatnonzero(np.abs(dif) > UMBRAL)
    return int(activo[0]), int(activo[-1])


def digitalizar(dif):
    """Diferencia de tension a bits logicos, con histeresis.

    En esta captura la diferencia positiva corresponde al nivel bajo.
    """
    estado = 0 if dif[0] > 0 else 1
    salida = np.zeros(len(dif), dtype=int)

    for i, v in enumerate(dif):
        if v > UMBRAL:
            estado = 0
        elif v < -UMBRAL:
            estado = 1
        salida[i] = estado

    return salida


def flancos(nivel, t):
    return t[np.flatnonzero(np.diff(nivel)) + 1]


def medir_bit(huecos):
    """El tiempo de bit, a partir de los huecos entre flancos.

    Con 8 muestras por bit cada flanco se localiza con un error de hasta una
    muestra, que es el 12% de un bit. Por eso el hueco mas corto no sirve:
    puede salir corto por puro redondeo. El mas *frecuente* si, porque el error
    se reparte a los dos lados y el valor central gana.
    """
    valores, cuentas = np.unique(np.round(huecos, 6), return_counts=True)
    bit = float(valores[np.argmax(cuentas)])

    for _ in range(6):                      # cada hueco dura un numero entero
        cuantos = np.round(huecos / bit)    # de bits
        cuantos[cuantos < 1] = 1
        encajan = np.abs(huecos / bit - cuantos) < 0.30
        bit = huecos[encajan].sum() / cuantos[encajan].sum()

    return float(bit)


def leer(t, nivel, instante):
    if instante < t[0] or instante > t[-1]:
        return None
    return int(nivel[int(np.searchsorted(t, instante))])


def buscar_flanco(momentos, esperado, margen):
    """El flanco real mas cercano al esperado, si cae dentro del margen."""
    if len(momentos) == 0:
        return None
    i = int(np.argmin(np.abs(momentos - esperado)))
    return momentos[i] if abs(momentos[i] - esperado) <= margen else None


def extraer_bytes(t, nivel, momentos, bit, inicio):
    """Recorre la trama byte a byte, resincronizando en cada BSS.

    Termina limpiamente al encontrar el FES: un bit bajo seguido de uno alto
    donde tocaria un BSS.
    """
    fin_tss = momentos[0]
    bits_tss = round((fin_tss - inicio) / bit)

    borde = fin_tss + bit                   # tras el TSS va el FSS, luego BSS
    octetos, avisos, cerrada = [], [], False

    while True:
        primero = leer(t, nivel, borde + 0.5 * bit)
        segundo = leer(t, nivel, borde + 1.5 * bit)

        if primero is None or segundo is None:
            avisos.append("la captura se acaba dentro de la trama")
            break

        if primero == 0 and segundo == 1:   # FES: fin de trama
            cerrada = True
            break

        if primero != 1 or segundo != 0:
            avisos.append(f"secuencia de byte rota en {borde:.2f} us")
            break

        bajada = buscar_flanco(momentos, borde + bit, 0.5 * bit)
        if bajada is None:
            avisos.append(f"flanco de sincronizacion perdido en {borde:.2f} us")
            break

        octeto = [leer(t, nivel, bajada + (k + 1.5) * bit) for k in range(8)]
        if None in octeto:
            avisos.append("la captura se acaba a mitad de un byte")
            break

        octetos.append(int("".join(str(b) for b in octeto), 2))
        borde = bajada + 9 * bit

    return bits_tss, octetos, cerrada, avisos


# ------------------------------------------------------------- protocolo

def a_bits(octetos):
    return [(b >> i) & 1 for b in octetos for i in range(7, -1, -1)]


def entero(bits):
    return int("".join(str(b) for b in bits), 2)


def crc(bits, polinomio, inicio, ancho):
    """Division binaria: el CRC es el resto."""
    registro = inicio
    mascara = (1 << ancho) - 1

    for b in bits:
        siguiente = b ^ ((registro >> (ancho - 1)) & 1)
        registro = (registro << 1) & mascara
        if siguiente:
            registro ^= polinomio

    return registro


def decodificar(octetos):
    """Interpreta los 40 bits de cabecera y verifica los dos CRC."""
    if len(octetos) < 8:
        return None

    cab = a_bits(octetos[:5])

    ident = entero(cab[5:16])
    palabras = entero(cab[16:23])
    crc_cab_leido = entero(cab[23:34])

    carga = octetos[5:5 + palabras * 2]
    cola = octetos[5 + palabras * 2:5 + palabras * 2 + 3]

    # CRC de cabecera: 11 bits sobre sincronismo, arranque, ID y longitud
    entrada_cab = cab[3:5] + cab[5:16] + cab[16:23]
    crc_cab_propio = crc(entrada_cab, CRC_CABECERA, INICIO_CABECERA, 11)

    # CRC de trama: 24 bits sobre la cabecera entera mas la carga util
    crc_trama_leido = entero(a_bits(cola)) if len(cola) == 3 else None
    entrada_trama = cab + a_bits(carga)

    canal = None
    crc_trama_propio = None
    for nombre, semilla in INICIO_CANAL.items():
        valor = crc(entrada_trama, CRC_TRAMA, semilla, 24)
        if valor == crc_trama_leido:
            canal, crc_trama_propio = nombre, valor
            break
    if crc_trama_propio is None:
        crc_trama_propio = crc(entrada_trama, CRC_TRAMA, INICIO_CANAL["A"], 24)

    return {
        "reservado": cab[0],
        "preambulo": cab[1],
        "no_nula": cab[2],
        "sincronismo": cab[3],
        "arranque": cab[4],
        "id": ident,
        "palabras": palabras,
        "bytes_carga": palabras * 2,
        "carga": carga,
        "ciclo": entero(cab[34:40]),
        "crc_cab_leido": crc_cab_leido,
        "crc_cab_propio": crc_cab_propio,
        "crc_cab_ok": crc_cab_leido == crc_cab_propio,
        "crc_trama_leido": crc_trama_leido,
        "crc_trama_propio": crc_trama_propio,
        "crc_trama_ok": crc_trama_leido == crc_trama_propio,
        "canal": canal,
    }


def analizar(ruta):
    t, dif = cargar(ruta)
    nivel = digitalizar(dif)

    a, b = acotar(dif)
    momentos = flancos(nivel[a:b + 1], t[a:b + 1])
    bit = medir_bit(np.diff(momentos))

    bits_tss, octetos, cerrada, avisos = extraer_bytes(t, nivel, momentos, bit, t[a])

    return {
        "muestras": len(t),
        "paso_us": float(t[1] - t[0]),
        "bit_us": bit,
        "mbits": 1.0 / bit,
        "muestras_por_bit": bit / (t[1] - t[0]),
        "bits_tss": bits_tss,
        "octetos": octetos,
        "cerrada": cerrada,
        "avisos": avisos,
        "trama": decodificar(octetos),
    }


# ------------------------------------------------------------- salida

def imprimir(r, ruta):
    print(f"\n{ruta}")
    print(f"  muestras            {r['muestras']}")
    print(f"  paso de muestreo    {r['paso_us']:.4f} us   ({1 / r['paso_us']:.0f} MS/s)")
    print(f"  tiempo de bit       {r['bit_us'] * 1000:.1f} ns")
    print(f"  velocidad del bus   {r['mbits']:.2f} Mbit/s")
    print(f"  muestras por bit    {r['muestras_por_bit']:.1f}")

    print(f"\n  TSS                 {r['bits_tss']} bits")
    print(f"  bytes leidos        {len(r['octetos'])}")
    print(f"  fin de trama        {'FES correcto' if r['cerrada'] else 'NO ENCONTRADO'}")

    for aviso in r["avisos"]:
        print(f"  aviso: {aviso}")

    m = r["trama"]
    if m is None:
        print("\n  trama incompleta: no se puede interpretar")
        return

    tipo = []
    if m["sincronismo"]:
        tipo.append("sincronismo")
    if m["arranque"]:
        tipo.append("arranque")
    if not m["no_nula"]:
        tipo.append("nula")

    print(f"\n  identificador       {m['id']}")
    print(f"  tipo                {', '.join(tipo) if tipo else 'datos'}")
    print(f"  contador de ciclo   {m['ciclo']}")
    print(f"  carga util          {m['palabras']} palabras = {m['bytes_carga']} bytes")

    for i in range(0, len(m["carga"]), 16):
        print(f"     {' '.join(f'{b:02X}' for b in m['carga'][i:i + 16])}")

    print(f"\n  CRC de cabecera     0x{m['crc_cab_leido']:03X} leido / "
          f"0x{m['crc_cab_propio']:03X} calculado   "
          f"{'CORRECTO' if m['crc_cab_ok'] else 'NO COINCIDE'}")
    print(f"  CRC de trama        0x{m['crc_trama_leido']:06X} leido / "
          f"0x{m['crc_trama_propio']:06X} calculado   "
          f"{'CORRECTO' if m['crc_trama_ok'] else 'NO COINCIDE'}")

    if m["canal"]:
        print(f"  canal               {m['canal']}   "
              f"(deducido de la semilla del CRC que cuadra)")


def graficar(ruta, salida="flexray.png"):
    import matplotlib.pyplot as plt

    tabla = pd.read_csv(ruta, skiprows=[1]).dropna()
    t = tabla.iloc[:, 0].to_numpy(dtype=float)
    bp = tabla.iloc[:, 1].to_numpy(dtype=float)
    bm = tabla.iloc[:, 2].to_numpy(dtype=float)

    fig, (arriba, abajo) = plt.subplots(2, 1, figsize=(13, 6), sharex=True)

    arriba.plot(t, bp, linewidth=0.7, label="BP")
    arriba.plot(t, bm, linewidth=0.7, label="BM")
    arriba.set_ylabel("V")
    arriba.legend(loc="upper right")
    arriba.set_title("Los dos hilos")

    abajo.plot(t, bp - bm, linewidth=0.7, color="black")
    abajo.axhline(UMBRAL, linestyle="--", linewidth=0.8)
    abajo.axhline(-UMBRAL, linestyle="--", linewidth=0.8)
    abajo.set_ylabel("V")
    abajo.set_xlabel("tiempo (us)")
    abajo.set_title("Diferencia BP menos BM, con los dos umbrales")

    for eje in (arriba, abajo):
        eje.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(salida, dpi=110)
    print(f"\n  grafica guardada en {salida}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", help="CSV exportado de PicoScope")
    parser.add_argument("--grafica", action="store_true")
    args = parser.parse_args()

    imprimir(analizar(args.csv), args.csv)

    if args.grafica:
        graficar(args.csv)


if __name__ == "__main__":
    main()
