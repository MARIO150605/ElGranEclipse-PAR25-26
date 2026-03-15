# Paradigmas de Programación Práctica I – Curso 2025/26

---

# El Gran Eclipse 

Como seguramente ya sabréis, el próximo 12 de agosto nuestra ciudad será testigo de un eclipse total de sol. Aunque estos fenómenos ocurren cada 18 meses de media, la estrechez de la trayectoria de la sombra lunar (entre 100 y 250 km) hace que desde un punto concreto solo se pueda ver un eclipse solar total cada 400 años (en promedio).

La excepcionalidad de este evento merece que dediquemos esta práctica al objetivo de crear una aplicación que nos permita conocer los eclipses que se van a producir en un determinado rango de fechas, los lugares del planeta desde los que se van a poder observar y pueda mostrar una animación de su desarrollo. 

Aunque en esta primera práctica la salida del programa va a realizarse exclusivamente en modo texto (consola), intentaremos exprimir al máximo las posibilidades de ese medio para que el resultado sea lo más informativo posible, mostrando para cada eclipse un mapa con su trayectoria en el que podamos hacer zoom en zonas concretas, y una animación en la que mostremos la evolución del eclipse cada 20 segundos. 

Para que este problema se pueda reducir a una complejidad adecuada vamos a apoyarnos en los siguientes elementos: 

- Usaremos la librería ephem para que se encargue de todos los cálculos astronómicos que nos permitan conocer con precisión la posición del sol y la luna en un momento y posición determinado. 

- Usaremos la librería global_land_mask para poder generar un mapa sencillo (solo distinguiremos entre mar y tierra) de una zona del planeta.

- Para los “gráficos” utilizaremos caracteres Unicode y códigos ANSI para asignarles color de fondo y primer plano. 
