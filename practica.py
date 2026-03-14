# *****************************************************************
#    PRACTICA 1 - ECLIPSES - PARADIGMAS DE PROGRAMACIÓN 2025-26
#    PLANTILLA DE LA PRÁCTICA
#    Apellidos, Nombre: 
#    Apellidos, Nombre: 
# *****************************************************************

# ephem: https://rhodesmill.org/pyephem/index.html
# ephem: https://pypi.org/project/ephem/
# ephem: pip install ephem

# global_land_mask: https://pypi.org/project/global-land-mask/
# global_land_mask: pip install global-land-mask
# global_land_mask: globe.is_land(lat, lon), lat-lon in degrees

# ansi: https://gist.github.com/fnky/458719343aabd01cfb17a3a4f7296797

import math
from ephem import *
from global_land_mask import globe
from time import sleep, time

# ************************ PARÁMETROS DEL PROGRAMA *****************************

N_SEG = 7                      # Número de segmentos en los que se divide cada intervalo para buscar el mínimo
PRECISION = (1.0/86400)/N_SEG  # Precisión (en fracciones de dia) en el cálculo de la fecha del eclipse (1 segundo)
SEP_LIM = 1.5 * degree         # Separación mínima para que pueda existir un eclipse de sol (radianes)
LIM_LAT = 70                   # Latitud máxima que se representa (grados)
N_FIL = 50                     # Número de filas de las matrices de mapa y animación
N_COL = 128                    # Número de columnas de las matrices de mapa y animación
ANIM_TAM_Y = 0.03              # Tamaño vertical (altura en radianes) del cuadro de animación
ANIM_TAM_X = ANIM_TAM_Y*N_COL/N_FIL  # Tamaño horizontal (azimut en radianes) del cuadro de animación

# COLORES ANSI-256 PARA MAPA (función col_mapa)

COLS_MAPA = ((18, 25, 32, 75),    # Gradiente agua
             (52, 88, 130, 178))  # Gradiente tierra
COL_SITIO = 196                   # Marca del sitio (rojo intenso)

# COLORES ANSI-256 PARA ANIMACIÓN (función col_anim)

COLS_ANIM = (75, 226, 242, 242, 130, 214, 172, 172)


# ************************ FUNCIONES PREDEFINIDAS *****************************

def cls():
    """ Borrado de pantalla (códigos ANSI) """
    print("\x1b[H\x1b[2J\x1b[3J")


def traduce_latlon(txt: str) -> (str | float, str | float):
    """ Traduce una posición geográfica en formato Google Maps (41°28'40.2"N 4°35'53.8"W) a una
        tupla de strings (latitud, longitud) con el formato de ephem
    :param txt: La posición obtenida de Google Maps
    :return: Una tupla de strings o floats (latitud, longitud)
    """
    if ',' in txt:
        return tuple(float(s) * degree for s in txt.split(", "))
    lat, lon = txt.split(' ')
    lat = lat.replace('°', ':').replace("'", ':')
    lon = lon.replace('°', ':').replace("'", ':')
    return (lat[0:-2] if lat[-1] == 'N' else '-' + lat[0:-2],
            lon[0:-2] if lon[-1] == 'E' else '-' + lon[0:-2])


def dist_ang(lon1: float, lat1: float, lon2: float, lat2: float) -> float:
    """  Calcula la distancia angular entre dos puntos en una esfera
    :param lon1: Longitud o Ascensión Recta o Azimut del primer punto (radianes)
    :param lat1: Latitud o Declinación o Altura del primer punto (radianes)
    :param lon2: Longitud o Ascensión Recta o Azimut del segundo punto (radianes)
    :param lat2: Latitud o Declinación o Altura del segundo punto (radianes)
    :return: Distancia angular (radianes) entre los puntos
    """
    return acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2))


def ocultacion(r1: float, r2: float, d: float) -> float:
    """ Calcula el porcentaje de area ocultada de un círculo por otro círculo
    :param r1: Radio del círculo ocultado
    :param r2: Radio del círculo ocultador
    :param d: Distancia entre los centros de ambos círculos
    :return: Un valor real entre 0 (no hay solapamiento) y 100 (totalmente ocultado)
    """
    if d > r1 + r2:  # Sin solapamiento
        return 0.0
    if d + r2 <= r1:  # Ocultador mas pequeño y completamente dentro del ocultado
        return 100 * r2 * r2 / (r1 * r1)
    if d + r1 <= r2:  # Ocultador mayor y tapa completamente al ocultado
        return 100.0
    # Area de la region de interseccion
    ai = r1 * r1 * acos((d * d + r1 * r1 - r2 * r2) / (2 * d * r1)) + \
        r2 * r2 * acos((d * d - r1 * r1 + r2 * r2) / (2 * d * r2)) - \
        0.5 * math.sqrt((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2))
    # Ratio del area del primer círculo
    return 100 * ai / (pi * r1 * r1)


def sep_alt_ocult(fec: Date, obj1: Body, obj2: Body, obs: Observer | None = None) -> (float, float, float):
    """ Calcula la distancia angular entre los cuerpos en esa fecha, la altura sobre el horizonte del primer cuerpo y
    el porcentaje de ocultación del primer cuerpo por parte del segundo.
    Si no se proporciona punto de observación los resultados son geocéntricos.
    :param fec: Fecha
    :param obj1: Cuerpo celeste ocultado
    :param obj2: Cuerpo celeste ocultador
    :param obs: Lugar de observación (None si el cálculo es geocéntrico)
    :return: Una tupla con la separación angular entre los cuerpos (radianes), altura sobre el horizonte (radianes) y
     porcentaje de ocultación (0-100)
    """
    if obs is None:
        # Cálculo geocéntrico, usamos RA/DEC y devolvemos una altura positiva cualquiera
        obj1.compute(fec)
        obj2.compute(fec)
        sep, alt = dist_ang(obj1.g_ra, obj1.g_dec, obj2.g_ra, obj2.g_dec), 0.00001
    else:
        # Cálculo topocéntrico, usamos Azimut/Altura
        obs.date = fec
        obj1.compute(obs)
        obj2.compute(obs)
        sep, alt = dist_ang(obj1.az, obj1.alt, obj2.az, obj2.alt), obj1.alt
    return sep, alt, ocultacion(obj1.radius, obj2.radius, sep)


def en_circulo(cx: float, cy: float, r: float, x: float, y: float):
    """ Comprueba si el punto (x,y) se encuentra en el círculo de centro (cx,cy) y radio r
    :param cx: Centro (x) del círculo
    :param cy: Centro (y) del círculo
    :param r: Radio del círculo
    :param x: Posición (x) del punto
    :param y: Posición (y) del punto
    :return: True si el punto está dentro del círculo
    """
    return (cx - x) * (cx - x) + (cy - y) * (cy - y) <= r * r


def fec_local(fec: Date) -> str:
    """ Convierte una fecha a texto en hora local """
    f = fec if isinstance(fec, Date) else Date(fec)
    return localtime(f).isoformat(' ', 'seconds') if f >= 25567.5 else str(f)


def interv(cen: float, tam: float, n: int) -> [float]:
    """ Devuelve una lista de los puntos centrales de un intervalo divido en segmentos contiguos iguales
    :param cen: Centro del intervalo
    :param tam: Longitud del intervalo
    :param n: Número de segmentos en que se divide el intervalo
    :return: Lista de posiciones de los puntos centrales de cada uno de los n segmentos, ordenada
    """
    x0 = cen - tam * (n - 1) / (2 * n)
    return [x0 + i * tam / n for i in range(n)]


def col_mapa(ocul: float, lat: float, lon: float) -> int:
    """ Devuelve el color de un pixel del mapa dado el nivel de ocultación y de si es tierra u océano
    :param ocul: El porcentaje (0-100) de ocultación
    :param lat: Latitud del punto (radianes)
    :param lon: Longitud del punto (radianes)
    :return: Color ANSI-256 del pixel
    """
    return COLS_MAPA[int(globe.is_land(lat/degree, lon/degree))][int(1.8 * math.log10(102 - ocul))]


def col_anim(en_obj1: bool, en_obj2: bool, bajo_horiz: bool = False) -> int:
    """ Devuelve el color del un pixel de animación según las condiciones del punto
    :param en_obj1: El punto está dentro del círculo que representa el primer cuerpo celeste
    :param en_obj2: El punto está dentro del círculo que representa el segundo cuerpo celeste
    :param bajo_horiz: El punto está bajo el horizonte
    :return: Color ANSI-256 del pixel
    """
    return COLS_ANIM[en_obj1 | en_obj2 << 1 | bajo_horiz << 2]


# ************************ FUNCIONES AUXILIARES DEL ALUMNO *******************************
def fecha_minima_precision(obj1, obj2, obs, fec_ini, duracion):
    """ Algoritmo 2.2.2: Refina la fecha del mínimo recursivamente """
    # 1. Dividir el intervalo en N_SEG segmentos (N_SEG+1 puntos)
    paso = duracion / N_SEG
    fechas = [Date(fec_ini + i * paso) for i in range(N_SEG + 1)]
    
    # Calcular separaciones y alturas
    datos = [sep_alt_ocult(f, obj1, obj2, obs) for f in fechas]
    
    # 2. Buscar el mínimo con altura >= 0
    # Si es geocéntrico (obs=None), alt siempre es > 0
    min_sep = float('inf')
    f_min_idx = -1
    
    for i, (sep, alt, ocu) in enumerate(datos):
        if alt >= 0 and sep < min_sep:
            min_sep = sep
            f_min_idx = i
            
    # Si todos están bajo el horizonte, devolver el centro del intervalo
    if f_min_idx == -1:
        return Date(fec_ini + duracion / 2)
    
    # 3. Si la duración es menor que la PRECISION, terminar
    if duracion < PRECISION:
        return fechas[f_min_idx]
    
    # 4. Recursividad: Nuevo intervalo alrededor del mínimo encontrado
    nueva_duracion = 2 * paso
    nueva_fec_ini = Date(fechas[f_min_idx] - paso)
    
    return fecha_minima_precision(obj1, obj2, obs, nueva_fec_ini, nueva_duracion)

def imprimir_tabla_eclipses(lista):
    """ Imprime la cabecera y las filas de los eclipses encontrados """
    print("\n  # | Fecha (Local)        |   Alt. |   Ocul. | Categoría")
    print("-" * 65)
    
    for i, e in enumerate(lista):
        # e[1] es la fecha topocéntrica (Local)
        # e[4] es la altura, e[5] la ocultación, e[6] la categoría
        print(f"  {i+1:d} | {fec_local(e[1])} | {e[4]:5.1f}º | {e[5]:6.2f}% | {e[6]}")
    
    print("-" * 65)
"""
def mostrar_mapa_eclipse(eclipse, obs):
    #Dibuja el mapa mundial con el rastro del eclipse para un eclipse seleccionado
    # Intervalos de coordenadas
    longitudes = interv(0, 2*math.pi, N_COL)           # [-180°,180°]
    latitudes = interv(0, -math.pi*LIM_LAT/90, N_FIL) # [+70° a -70°]

    # Recorrer de 2 en 2 filas (para medio bloque)
    for i in range(0, N_FIL, 2):
        lat_sup = latitudes[i]
        lat_inf = latitudes[i+1] if i+1 < N_FIL else latitudes[i]

        for j in range(N_COL):
            lon = longitudes[j]

            # Fila superior
            obs.lat, obs.lon = lat_sup, lon
            fecha_ini = Date(eclipse[0]-1)  # Día anterior a fecha geocéntrica
            f_min_sup = fecha_minima_precision(Sun(), Moon(), obs, fecha_ini, 1.0)
            _, _, ocul_sup = sep_alt_ocult(f_min_sup, Sun(), Moon(), obs)
            color_sup = col_mapa(ocul_sup, lat_sup, lon)

            # Fila inferior
            obs.lat, obs.lon = lat_inf, lon
            f_min_inf = fecha_minima_precision(Sun(), Moon(), obs, fecha_ini, 1.0)
            _, _, ocul_inf = sep_alt_ocult(f_min_inf, Sun(), Moon(), obs)
            color_inf = col_mapa(ocul_inf, lat_inf, lon)

            # Imprimir medio bloque
            print(f"\x1b[38;5;{color_sup}m\x1b[48;5;{color_inf}m▀", end="")

        print("\x1b[0m")  # reset colores al final de la línea
"""
def mostrar_mapa_eclipse_zoom(eclipse, obs: Observer, lon_cen, lat_cen, lon_tam, lat_tam):
    """Similar a mostrar_mapa_eclipse, pero permite zoom en cualquier zona"""
    longitudes = interv(lon_cen, lon_tam, N_COL)
    latitudes = interv(lat_cen, -lat_tam, N_FIL)

    sitio_lat = obs.lat
    sitio_lon = obs.lon
    
    res_lon = lon_tam / N_COL
    res_lat = lat_tam / N_FIL
    
    for i in range(0, N_FIL, 2):
        lat_sup = latitudes[i]
        lat_inf = latitudes[i+1] if i+1 < N_FIL else latitudes[i]

        for j in range(N_COL):
            lon = longitudes[j]

            # fila superior
            obs.lat, obs.lon = lat_sup, lon
            fecha_ini = Date(eclipse[0]-1)
            f_min_sup = fecha_minima_precision(Sun(), Moon(), obs, fecha_ini, 1.0)
            _, _, ocul_sup = sep_alt_ocult(f_min_sup, Sun(), Moon(), obs)
            
            # ¿Es este píxel la ubicación del sitio?
            if abs(lat_sup - sitio_lat) < res_lat/2 and abs(lon - sitio_lon) < res_lon/2:
                color_sup = COL_SITIO
            else:
                color_sup = col_mapa(ocul_sup, lat_sup, lon)
                
            #color_sup = col_mapa(ocul_sup, lat_sup, lon)

            # fila inferior
            obs.lat, obs.lon = lat_inf, lon
            f_min_inf = fecha_minima_precision(Sun(), Moon(), obs, fecha_ini, 1.0)
            _, _, ocul_inf = sep_alt_ocult(f_min_inf, Sun(), Moon(), obs)
            
            if abs(lat_inf - sitio_lat) < res_lat/2 and abs(lon - sitio_lon) < res_lon/2:
                color_inf = COL_SITIO
            else:
                color_inf = col_mapa(ocul_inf, lat_inf, lon)
                
            #color_inf = col_mapa(ocul_inf, lat_inf, lon)

            # imprimir medio bloque
            print(f"\x1b[38;5;{color_sup}m\x1b[48;5;{color_inf}m▀", end="")
        print("\x1b[0m")
    
    obs.lat = sitio_lat
    obs.lon = sitio_lon
    # Restaurar la posición original del observador (modificada durante el cálculo del mapa)

def animacion_eclipse(eclipse, obs: Observer, obj1: Body, obj2: Body):
    """ Muestra una animación del eclipse desde el punto de observación
    :param eclipse: Tupla con los datos del eclipse
    :param obs: Punto de observación
    :param obj1: Cuerpo celeste ocultado (Sol)
    :param obj2: Cuerpo celeste ocultador (Luna)
    """
    # 1. Buscar el inicio: partimos del momento de máxima ocultación 
    #    y restamos 1 segundo hasta que la ocultación sea 0
    fec_act = Date(eclipse[1])
    _, _, ocu = sep_alt_ocult(fec_act, obj1, obj2, obs)
    while ocu > 0:
        fec_act = Date(fec_act - 1.0 / 86400)
        _, _, ocu = sep_alt_ocult(fec_act, obj1, obj2, obs)

    # 2. Borrar pantalla y mostrar el primer fotograma 
    cls()
    dibujar_fotograma(fec_act, ocu, obs, obj1, obj2)

    # 3. Bucle: avanzar 20s, calcular ocultación, mostrar fotograma, parar cuando ocu=0
    paso_20s = 20.0 / 86400
    while True:
        fec_act = Date(fec_act + paso_20s)
        _, _, ocu = sep_alt_ocult(fec_act, obj1, obj2, obs)
        dibujar_fotograma(fec_act, ocu, obs, obj1, obj2)
        if ocu <= 0:
            break


def dibujar_fotograma(fec_act, ocu, obs, obj1, obj2):
    """ Dibuja un fotograma de la animación del eclipse
    :param fec_act: Fecha actual del fotograma
    :param ocu: Porcentaje de ocultación en este fotograma
    :param obs: Punto de observación
    :param obj1: Sol (ya computado por sep_alt_ocult)
    :param obj2: Luna (ya computado por sep_alt_ocult)
    """
    # Listas de alturas (mayor a menor) y azimuts centradas en el Sol
    l_alt = interv(0, -ANIM_TAM_Y, N_FIL)
    l_azi = interv(0,  ANIM_TAM_X, N_COL)

    # Posición absoluta del Sol y la Luna (actualizados por sep_alt_ocult)
    a_sol = float(obj1.alt)
    z_sol = float(obj1.az)
    a_lun = float(obj2.alt)
    z_lun = float(obj2.az)
    r_sol = float(obj1.radius)
    r_lun = float(obj2.radius)

    # Posición relativa de la Luna respecto al Sol (Sol en (0,0))
    rel_l_az  = z_lun - z_sol
    rel_l_alt = a_lun - a_sol
    # Horizonte relativo: altitud 0 absoluta = -a_sol relativo al Sol
    rel_h_alt = -a_sol

    # Construir el buffer: cursor a (0,0), fecha y ocultación, luego la matriz
    buf = f"\x1b[HFecha: {fec_local(fec_act)} | Ocultación: {ocu:6.2f}%\x1b[K\n"

    # Compresión de filas: dos filas de la matriz en una línea de texto con ▀
    for i in range(0, N_FIL, 2):
        fila = ""
        y_s = l_alt[i]
        y_i = l_alt[i + 1] if i + 1 < N_FIL else y_s

        for x in l_azi:
            # Semipíxel superior
            es_s = en_circulo(0, 0, r_sol, x, y_s)
            el_s = en_circulo(rel_l_az, rel_l_alt, r_lun, x, y_s)
            bh_s = y_s < rel_h_alt
            # Semipíxel inferior
            es_i = en_circulo(0, 0, r_sol, x, y_i)
            el_i = en_circulo(rel_l_az, rel_l_alt, r_lun, x, y_i)
            bh_i = y_i < rel_h_alt

            fila += f"\x1b[38;5;{col_anim(es_s, el_s, bh_s)}m\x1b[48;5;{col_anim(es_i, el_i, bh_i)}m▀"

        buf += fila + "\x1b[0m\x1b[K\n"

    print(buf, end="", flush=True)
    sleep(0.02)


# ************************ CABECERAS DE FUNCIONES AUXILIARES *****************************

def lista_eclipses(tipo: int, fec_ini: Date, fec_fin: Date, obj1: Body, obj2: Body, obs: Observer, min_ocu: float)\
        -> [(Date, Date, float, float, float, float, str)]:
    """ Devuelve la lista de los eclipses que se pueden dar entre dos fechas, filtrando por separación
    mínima desde el punto de vista geocéntrico y ocultación mínima desde el punto de vista del observador
    :param tipo: Tipo de eclipse (1 -> Eclipse de Sol)
    :param fec_ini: Fecha inicial
    :param fec_fin: Fecha final
    :param obj1: Cuerpo celeste ocultado
    :param obj2: Cuerpo celeste ocultador
    :param obs: Punto de observación
    :param min_ocu: Mínimo valor de ocultación para ser incluido (porcentaje, 1-100)
    :return: Enumeración de datos de eclipses: (fecha geo, fecha top, sep geo, sep top, altitud, ocultación, mensaje),
     las separaciones y altitud en grados.
    """
    res = []
    # Empezamos con tres fechas consecutivas (paso de 1 día)
    f1 = Date(fec_ini)
    f2 = Date(f1 + 1)
    f3 = Date(f2 + 1)
    
    while f3 <= fec_fin:
        # 1. Calcular separaciones geocéntricas (s1, s2, s3)
        s1, _, _ = sep_alt_ocult(f1, obj1, obj2, None)
        s2, _, _ = sep_alt_ocult(f2, obj1, obj2, None)
        s3, _, _ = sep_alt_ocult(f3, obj1, obj2, None)
        
        # 2. ¿Hay un mínimo local en f2? (s1 > s2 y s2 < s3)
        if s1 > s2 and s2 < s3:
            # 3. Fecha precisa geocéntrica (Intervalo de 48h desde f1)
            f_geo_min = fecha_minima_precision(obj1, obj2, None, f1, 2.0)
            s_g_min, _, _ = sep_alt_ocult(f_geo_min, obj1, obj2, None)
            
            # 4. Si la separación es <= SEP_LIM (1.5º), es un eclipse real
            if s_g_min <= SEP_LIM:
                # 5. Fecha precisa topocéntrica
                f_topo_min = fecha_minima_precision(obj1, obj2, obs, f1, 2.0)
                sep_t, alt_t, ocu_t = sep_alt_ocult(f_topo_min, obj1, obj2, obs)
                
                # 6. Filtrar por ocultación mínima
                if ocu_t >= min_ocu:
                    if ocu_t < 1.0: cat = "NO VISIBLE"
                    elif ocu_t >= 95.0: cat = "TOTAL"
                    else: cat = "PARCIAL"
                    
                    res.append((f_geo_min, f_topo_min, s_g_min/degree, 
                                sep_t/degree, alt_t/degree, ocu_t, cat))
            
            # Tras encontrarlo, avanzamos para no repetir el mismo eclipse
            f1 = Date(f3 + 25) # Saltamos casi un mes
        else:
            f1 = Date(f1 + 1) # Avanzamos un día (Punto 7 del PDF)
            
        f2 = Date(f1 + 1)
        f3 = Date(f2 + 1)
            
    return res


def menu_principal(tipo: int, lis: [(Date, Date, float, float, float, float, str)],
                   obj1: Body, obj2: Body, obs: Observer):
    """ Muestra el menu de acciones (listado, mapa y animación) sobre los eclipses
    :param tipo: Tipo de cálculo (0 -> Eclipses de Sol)
    :param lis: Lista de posibles eclipses y su información asociada
    :param obj1: Cuerpo celeste ocultado
    :param obj2: Cuerpo celeste ocultador
    :param obs: Lugar de observación
    """
    """ Muestra el menú de acciones sobre los eclipses encontrados """
    
    while True:
        # Mostrar la lista de eclipses
        imprimir_tabla_eclipses(lis)
        
        # Elegir eclipse
        seleccion = input(f"Elija un eclipse (1-{len(lis)}) o (S) para salir: ").strip().lower()
        if seleccion == 's':
            break
        
        if not seleccion.isdigit() or not (1 <= int(seleccion) <= len(lis)):
            print("Opción inválida, inténtelo de nuevo.")
            continue
        
        idx = int(seleccion) - 1
        eclipse = lis[idx]

        # Inicializar zoom
        lon_cen, lat_cen = 0, 0
        lon_tam, lat_tam = 2*math.pi, math.pi*LIM_LAT/90
        zoom_hist = []

        # Menú de visualización del mapa
        while True:
            # Mostrar mapa con zoom actual
            mostrar_mapa_eclipse_zoom(eclipse, obs, lon_cen, lat_cen, lon_tam, lat_tam)

            # Pedir acción al usuario
            cmd = input("\nOpciones: Zoom(1-9), Deshacer(0), Ant/Sig Eclipse(< >), Animacion(A), Volver(V): ").strip().lower()
            
            if cmd == 'v':  # Volver a la lista
                break
            elif cmd == '0':  # Deshacer zoom
                if zoom_hist:
                    lon_cen, lat_cen, lon_tam, lat_tam = zoom_hist.pop()
            elif cmd in '123456789':  # Zoom a una zona
                zona = int(cmd) - 1
                zoom_hist.append((lon_cen, lat_cen, lon_tam, lat_tam))
                lon_cen += (zona % 3 - 1)*(lon_tam / 4)
                lat_cen -= (zona // 3 - 1)*(lat_tam / 4) # esto se debe poner un - ya que las latitudes estan invertidas
                lon_tam /= 2
                lat_tam /= 2
            elif cmd == '<':  # Eclipse anterior
                if idx > 0:
                    idx -= 1
                    eclipse = lis[idx]
            elif cmd == '>':  # Eclipse siguiente
                if idx < len(lis)-1:
                    idx += 1
                    eclipse = lis[idx]
            elif cmd == 'a':  # Animación (pendiente)
                animacion_eclipse(eclipse,obs,obj1,obj2)
                input("Pulse retorno para continuar. ")
                continue
            else:
                print("Opción inválida, inténtelo de nuevo.")

# FUNCIÓN PRINCIPAL

def main():
    print("ECLIPSES - PARADIGMAS 2025-26")
    print("1. Eclipse de Sol")
    tipo = input("Escoja tipo de cálculo (1): ")
    tipo = int(tipo) if tipo else 1

    # Datos necesarios
    fec_ini = input("Fecha inicial (formato año-mes-dia): ")
    fec_fin = input("Fecha final (formato año-mes-dia):   ")
    loc_geo = input("Localidad, posicion geografica: ")
    loc_alt = input("Localidad, altura (m): ")
    min_ocu = input("Porcentaje minimo de ocultacion (0..100): ")
    # Traducción y sustitución de valores por defecto
    fec_ini = Date(fec_ini) if fec_ini else now()
    fec_fin = Date(fec_fin) if fec_fin else Date(fec_ini+3653)
    # La localidad por defecto es Valladolid
    lat, lon = traduce_latlon(loc_geo) if loc_geo else (41.66308134*degree, -4.70494676*degree)
    loc_alt = int(loc_alt) if loc_alt else 0
    min_ocu = int(min_ocu) if min_ocu else 0
    # Creación del punto de observación y de los cuerpos celestes
    sitio = Observer()
    sitio.lat = lat
    sitio.lon = lon
    sitio.elevation = loc_alt
    obj1, obj2 = Sun(), Moon()
    # Obtención de la lista de eclipses
    lis = lista_eclipses(tipo, fec_ini, fec_fin, obj1, obj2, sitio, min_ocu)
    # Menu principal del programa
    menu_principal(tipo, lis, obj1, obj2, sitio)


if __name__ == '__main__':
    main()
