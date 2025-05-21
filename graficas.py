from pylab import *
from matplotlib import style
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

def grafica3d(x,y,z,t,name,method):
    """Función para graficar los atractores en el espacio de fases y guardar en pdf"""
    plt.figure()
    plt.plot(x, y, "m", lw=0.3)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.savefig(name + "_" + method + "_xy_atractor.pdf", dpi=300, bbox_inches='tight')
    plt.clf()

    plt.figure()
    plt.plot(x, z, "m", lw=0.3)
    plt.xlabel("x")
    plt.ylabel("z")
    plt.savefig(name + "_" + method + "_xz_atractor.pdf", dpi=300, bbox_inches='tight')
    plt.clf()

    plt.figure()
    plt.plot(y, z, "m", lw=0.3)
    plt.xlabel("y")
    plt.ylabel("z")
    plt.savefig(name + "_" + method + "_yz_atractor.pdf", dpi=300, bbox_inches='tight')
    plt.clf()

    fig = plt.figure()
    ax1 = fig.add_subplot(111,projection='3d')
    ax1.plot(x,y,z,"m",lw=0.2)
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_zlabel('z')
    ax1.view_init(elev=-30, azim=-45) # Elevación y azimut en grados
    plt.savefig(name + "_" + method + "_" + "atractor3d_" + ".pdf",dpi=300,bbox_inches= 'tight')
    #plt.show()
    clf()
    
def grafica4d(x,y,z,w,t,name,method):
    """Función para graficar los atractores en el espacio de fases y guardar en pdf"""
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

    plt.savefig(name + "_atractores_" + ".pdf",dpi=300,bbox_inches= 'tight')
    clf()