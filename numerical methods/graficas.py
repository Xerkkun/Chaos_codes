from pylab import *
from matplotlib import style
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

def grafica3d(x,y,z,t,name,method):
    """Función para graficar los atractores en el espacio de fases y guardar en png"""
    plt.figure()
    plt.plot(x, y, "m", lw=0.3)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.savefig(name + "_" + method + "_xy_atractor.png", dpi=300, bbox_inches='tight')
    plt.clf()

    plt.figure()
    plt.plot(x, z, "m", lw=0.3)
    plt.xlabel("x")
    plt.ylabel("z")
    plt.savefig(name + "_" + method + "_xz_atractor.png", dpi=300, bbox_inches='tight')
    plt.clf()

    plt.figure()
    plt.plot(y, z, "m", lw=0.5)
    plt.xlabel("y")
    plt.ylabel("z")
    plt.savefig(name + "_" + method + "_yz_atractor.png", dpi=300, bbox_inches='tight')
    plt.clf()
    
    plt.figure(figsize=(10,5))
    plt.plot(t, x, "m", lw=0.5)
    plt.xlabel("t")
    plt.ylabel("x")
    plt.savefig(name + "_" + method + "_x.png", dpi=300, bbox_inches='tight')
    plt.clf()

    plt.figure(figsize=(10,5))
    plt.plot(t, y, "m", lw=0.5)
    plt.xlabel("t")
    plt.ylabel("y")
    plt.savefig(name + "_" + method + "_y.png", dpi=300, bbox_inches='tight')
    plt.clf()

    plt.figure(figsize=(10,5))
    plt.plot(t, z, "m", lw=0.3)
    plt.xlabel("t")
    plt.ylabel("z")
    plt.savefig(name + "_" + method + "_z.png", dpi=300, bbox_inches='tight')
    plt.clf()

    fig = plt.figure()
    ax1 = fig.add_subplot(111,projection='3d')
    ax1.plot(x,y,z,"m",lw=0.2)
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_zlabel('z')
    ax1.view_init(elev=-30, azim=-45) # Elevación y azimut en grados
    plt.savefig(name + "_" + method + "_" + "atractor3d_" + ".png",dpi=300,bbox_inches= 'tight')
    #plt.show()
    clf()
    
def grafica4d(x,y,z,w,t,name,method):
    """Función para graficar los atractores en el espacio de fases y guardar en png"""
    subplot(221)
    p1,=plot(x,y,"m",lw=0.3)
    xlabel("x")
    ylabel("y")
    #plt.axis('square')

    subplot(222)
    p2,=plot(y,z,"m",lw=0.3)
    xlabel("y")
    ylabel("z")

    subplot(223)
    p3,=plot(x,z,"m",lw=0.3)
    xlabel("x")
    ylabel("z")

    subplot(224)
    p4,=plot(x,w,"m",lw=0.3)
    xlabel("x")
    ylabel("w")

    plt.savefig(name + "_" + method +  "_atractores_" + ".png",dpi=300,bbox_inches= 'tight')
    clf()
    
def graph_colors(x,y,z,t,name,method):
    
    # Normalize the x-axis to use for color mapping
    norm = plt.Normalize(x.min(), x.max())
    colors = plt.cm.hsv(norm(x))

    # Set up the figure and axis with black background
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_facecolor('black')

    # Plot the Lorenz attractor with a color gradient dependent on the x position
    # Here we plot only the XY projection by setting the Z values to 0
    sc = ax.scatter(x, y, c=colors, s=0.5, alpha=0.6)

    # Hide the axes
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')  # Turn off the axis completely

    # Tight layout and save to png with a black background
    plt.tight_layout()
    plt.savefig(name + "_" + method +  "_xy_atractor.png", dpi=300, format='png', bbox_inches='tight', facecolor='black', edgecolor='none')

    # Close the plot
    plt.close()
#==================================================================
    # Normalize the x-axis to use for color mapping
    norm = plt.Normalize(x.min(), x.max())
    colors = plt.cm.hsv(norm(x))

    # Set up the figure and axis with black background
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_facecolor('black')

    # Plot the Lorenz attractor with a color gradient dependent on the x position
    # Here we plot only the XY projection by setting the Z values to 0
    sc = ax.scatter(x, z, c=colors, s=0.5, alpha=0.6)

    # Hide the axes
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')  # Turn off the axis completely

    # Tight layout and save to png with a black background
    plt.tight_layout()
    plt.savefig(name + "_" + method +  "_xz_atractor.png", dpi=300, format='png', bbox_inches='tight', facecolor='black', edgecolor='none')

    # Close the plot
    plt.close()
#==================================================================
      # Normalize the y-axis to use for color mapping
    norm = plt.Normalize(y.min(), y.max())
    colors = plt.cm.hsv(norm(y))

    # Set up the figure and axis with black background
    fig, ay = plt.subplots(figsize=(8, 8))
    ay.set_facecolor('black')

    # Plot the Lorenz attractor with a color gradient dependent on the x position
    # Here we plot only the XY projection by setting the Z values to 0
    sc = ay.scatter(y, z, c=colors, s=0.5, alpha=0.6)

    # Hide the axes
    ay.set_xticks([])
    ay.set_yticks([])
    ay.axis('off')  # Turn off the axis completely

    # Tight layout and save to png with a black background
    plt.tight_layout()
    plt.savefig(name + "_" + method +  "_yz_atractor.png", dpi=300, format='png', bbox_inches='tight', facecolor='black', edgecolor='none')

    # Close the plot
    plt.close()
    
def graph_colors_3d(x, y, z, t, name, method):
    # Normalizar el eje x para usarlo en el mapeo de colores
    norm = plt.Normalize(x.min(), x.max())
    colors = plt.cm.hsv(norm(x))

    # Configurar la figura y el eje con fondo negro
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_facecolor('black')

    # Graficar el atractor con un gradiente de color dependiente de la posición x
    sc = ax.scatter(x, y, z, c=colors, s=0.5, alpha=0.6)

    # Ocultar los ejes
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.axis('off')  # Desactivar completamente el eje
    ax.view_init(elev=-30, azim=-45)  # Elevación y azimut en grados
    # Ajustar el layout y guardar en PDF con fondo negro
    plt.tight_layout()
    plt.savefig(name + "_" + method + "_xyz_atractor.png", dpi=300, format='png', bbox_inches='tight', facecolor='black', edgecolor='none')

    # Cerrar la gráfica
    plt.close() 

    
def graph_optimization(sumax, sumay, sumaz, num_peaksx, num_peaksy, num_peaksz, num_peaks_thrx, num_peaks_thry, num_peaks_thrz, eigenvalues, DKY_mean, t_step, X, sol, t2, sol2, Y2, Y3, Y4, frq2 ):
    
    #Graficar 
    label_potx = 'pot=' + str(np.round(sumax,2)) + ',peaks=' + str(num_peaksx) + ',' + str(num_peaks_thrx)
    label_poty = 'pot=' + str(np.round(sumay,2)) + ',peaks=' + str(num_peaksy) + ',' + str(num_peaks_thry)
    label_potz = 'pot=' + str(np.round(sumaz,2)) + ',peaks=' + str(num_peaksz) + ',' + str(num_peaks_thrz)
    label_eig = str(np.round(eigenvalues.real, 2) + np.round(eigenvalues.imag, 2) * 1j)
    label_dky = str(np.round(-DKY_mean,4))
    label_h = str(t_step)
    title = str(np.round(X,2))
    #file_name = str(j) + ".png"
    
    fig = plt.figure(figsize=(10, 8))

    ax1 = fig.add_subplot(321)
    ax1.plot(sol[:,0], sol[:,2])
    plt.legend([label_dky])
    plt.xlabel('Eje x')
    plt.ylabel('Eje z')
    plt.title(title)
    plt.grid(True)

    ax2 = fig.add_subplot(322)
    ax2.scatter(eigenvalues.real, eigenvalues.imag, marker="o")
    plt.legend([label_eig])
    plt.xlabel('Re')
    plt.ylabel('Im')
    plt.axhline(0, color="black")
    plt.axvline(0, color="black")

    ax1 = fig.add_subplot(323)
    ax1.plot(t2, sol2[:,0])
    plt.legend([label_h])
    ax1.set_xlabel('Tiempo (s)')
    ax1.set_ylabel('$x(t)$')

    ax2 = fig.add_subplot(324)
    ax2.vlines(frq2, 0, np.abs(Y2.imag))
    plt.legend([label_potx])
    plt.xlim(-10, 10)
    plt.xlabel('Frecuencia (Hz)')
    plt.ylabel('Im($Y_x$)')
    
    ax2 = fig.add_subplot(325)
    ax2.vlines(frq2, 0, np.abs(Y3.imag))
    plt.legend([label_poty])
    plt.xlim(-10, 10)
    plt.xlabel('Frecuencia (Hz)')
    plt.ylabel('Im($Y_y$)')
    
    ax2 = fig.add_subplot(326)
    ax2.vlines(frq2, 0, np.abs(Y4.imag))
    plt.legend([label_potz])
    plt.xlim(-10, 10)
    plt.xlabel('Frecuencia (Hz)')
    plt.ylabel('Im($Y_z$)')
    
    fig.tight_layout()
    plt.show()

def temporal_evolution(x, y, z, t, t_step, transient, name, method):
    label_h = r'$h =$'+str(t_step)
    plt.figure(figsize=(30, 12))

    plt.plot(t, x, "m", lw=1.5 )
    plt.xlabel(r'$t$', fontsize=24)
    plt.ylabel(r'$x(t)$', fontsize=24)
    plt.legend([label_h], fontsize=24)
    # Cambiar el tamaño de los números de los ejes
    # Cambia '12' por el tamaño deseado para el eje x
    plt.tick_params(axis='x', labelsize=24)
    # Cambia '12' por el tamaño deseado para el eje y
    plt.tick_params(axis='y', labelsize=24)

    # Trazar una línea vertical en x=1
    # 'k' es el color negro, '--' es el estilo de línea discontinua, lw es el grosor de la línea
    plt.axvline(x=transient, color='k', linestyle='--', lw=2)
    
    # Ajusta esto para mover el texto más a la derecha
    # x_pos_texto = 1 + (plt.xlim()[1] - plt.xlim()[0]) * 0.65
    x_pos_texto = 1 + (plt.xlim()[1] - plt.xlim()[0]) * 0.03

    # Ajusta la posición vertical del texto
    y_pos_texto = plt.ylim()[0] + (plt.ylim()[1] - plt.ylim()[0]) * 0.05
    # plt.text(x_pos_texto, y_pos_texto, r'$t = 2\lambda_{max}$', verticalalignment='bottom',
    #          horizontalalignment='left', fontsize=20, color='k')


    # Guardar la figura
    plt.savefig(name + "_" + method + "_temporal_" +
                ".pdf", dpi=300, bbox_inches='tight')

    # Limpiar la figura para evitar superposiciones si se vuelve a llamar a plt.show() o se guarda otra figura
    plt.clf()


    clf()
