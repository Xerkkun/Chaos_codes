with open("salida_pso.lyaps", "r") as file:
    file.seek(0,2)
    file.seek(file.tell()-9)
    DKY = float(file.read())
    file.seek(0,2)
    file.seek(file.tell()-309)
    LE = float(file.read(12))

print(DKY,LE)