"""FlexRay: de la senal del osciloscopio a los bytes de la trama.

FlexRay no usa relleno de bits como CAN. Usa otra cosa: delante de CADA byte
manda dos bits de servicio (uno alto y uno bajo) llamados BSS. El flanco de
bajada del BSS es un punto de sincronizacion, y el receptor lo aprovecha para
recolocarse cada 10 bits. Aqui se hace lo mismo.

    python flexray_bits.py FlexRay_Trace.csv
"""

import sys

import numpy as np
import pandas as pd

UMBRAL = 0.30      # V en la diferencia entre los dos hilos; en medio se conserva


def cargar(ruta):
    tabla = pd.read_csv(ruta, skiprows=[1]).dropna()
    t = tabla.iloc[:, 0].to_numpy(dtype=float)
    bp = tabla.iloc[:, 1].to_numpy(dtype=float)
    bm = tabla.iloc[:, 2].to_numpy(dtype=float)
    return t, bp - bm


def acotar(dif):
    """Donde empieza y acaba la trama.

    En reposo los dos hilos estan al mismo potencial y la diferencia es casi
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
    """Instantes de cada cambio de nivel."""
    return t[np.flatnonzero(np.diff(nivel)) + 1]


def medir_bit(huecos):
    """El tiempo de bit, a partir de los huecos entre flancos.

    Con solo 8 muestras por bit, cada flanco se detecta con un error de hasta
    una muestra, que aqui es el 12% de un bit. Por eso el hueco mas corto no
    sirve: puede salir corto por puro redondeo. El hueco mas *frecuente* si,
    porque el error se reparte a los dos lados y el valor central gana.
    """
    valores, cuentas = np.unique(np.round(huecos, 6), return_counts=True)
    bit = float(valores[np.argmax(cuentas)])

    # Afinado: cada hueco dura un numero entero de bits.
    for _ in range(6):
        cuantos = np.round(huecos / bit)
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


def extraer(t, nivel, momentos, bit, inicio_tss):
    """Recorre la trama byte a byte, resincronizando en cada BSS."""
    # El TSS es una tirada de nivel bajo. Acaba en el primer flanco.
    fin_tss = momentos[0]
    bits_tss = round((fin_tss - inicio_tss) / bit)

    # Tras el TSS viene el FSS (un bit alto) y ya empieza el primer BSS.
    borde = fin_tss + bit
    bytes_leidos, avisos = [], []

    while True:
        # BSS: un bit alto y uno bajo. El flanco de bajada esta un bit despues.
        bajada = buscar_flanco(momentos, borde + bit, 0.5 * bit)

        if bajada is None:
            if leer(t, nivel, borde + 0.5 * bit) == 1:
                break                      # ya no hay BSS: fin de trama
            avisos.append(f"BSS perdido en {borde:.2f} us")
            break

        if leer(t, nivel, borde + 0.5 * bit) != 1:
            avisos.append(f"BSS con nivel raro en {borde:.2f} us")
            break

        # Los ocho bits del byte empiezan en el flanco de bajada del BSS.
        octeto = [leer(t, nivel, bajada + (k + 1.5) * bit) for k in range(8)]
        if None in octeto:
            break

        bytes_leidos.append(int("".join(str(b) for b in octeto), 2))
        borde = bajada + 9 * bit           # el BSS del byte siguiente

    return bits_tss, bytes_leidos, avisos


def main():
    ruta = sys.argv[1] if len(sys.argv) > 1 else "FlexRay_Trace.csv"

    t, dif = cargar(ruta)
    nivel = digitalizar(dif)

    # Solo interesa el tramo con trama: fuera de el la senal esta en reposo
    # y cualquier flanco que se detecte ahi es ruido.
    a, b = acotar(dif)
    momentos = flancos(nivel[a:b + 1], t[a:b + 1])
    bit = medir_bit(np.diff(momentos))

    print(f"\n{ruta}")
    print(f"  muestras            {len(t)}")
    print(f"  paso de muestreo    {t[1] - t[0]:.4f} us   ({1 / (t[1] - t[0]):.0f} MS/s)")
    print(f"  tiempo de bit       {bit * 1000:.1f} ns")
    print(f"  velocidad del bus   {1 / bit:.2f} Mbit/s")
    print(f"  muestras por bit    {bit / (t[1] - t[0]):.1f}")

    bits_tss, octetos, avisos = extraer(t, nivel, momentos, bit, t[a])

    print(f"\n  TSS                 {bits_tss} bits de nivel bajo")
    print(f"  bytes leidos        {len(octetos)}")

    for aviso in avisos:
        print(f"  aviso: {aviso}")

    if not octetos:
        return

    print("\n  cabecera (5 bytes)")
    print(f"     {' '.join(f'{b:02X}' for b in octetos[:5])}")

    if len(octetos) > 8:
        carga = octetos[5:-3]
        print(f"\n  carga util ({len(carga)} bytes)")
        for i in range(0, len(carga), 16):
            print(f"     {' '.join(f'{b:02X}' for b in carga[i:i + 16])}")

        print("\n  cola: CRC de trama (3 bytes)")
        print(f"     {' '.join(f'{b:02X}' for b in octetos[-3:])}")


if __name__ == "__main__":
    main()
