import vtk

def write_points(points,
                 fh,
                 index):
    for idx, i in enumerate(range(0, points.GetNumberOfPoints())):
        x = points.GetPoint(i)
        fh.write(str(x[index]))
        if idx % 5 == 0 and idx != 0:
            fh.write("\n")
        else:
            fh.write(" ")
    fh.write("\n")

def write_array(array,
                fh):
    for idx, i in enumerate(range(0, array.GetNumberOfTuples())):
        fh.write(str(array.GetTuple1(i)))
        if idx % 5 == 0 and idx != 0:
            fh.write("\n")
        else:
            fh.write(" ")
    fh.write("\n")

def write_tecplot_ascii(ug,
                        fh,
                        title,
                        time):
                         
    if not ug.IsHomogeneous():
        print "WARNING: Unable to convert non-homogenous (cells) VTK unstructured grid"
        return
 
    if not ug.GetNumberOfCells():
        return

    ct = ug.GetCell(0).GetCellType()
    is3d = True
    if ct == vtk.VTK_LINE: 
        tecplot_type = "FELINESEG"
        is3d = False
    elif ct == vtk.VTK_TRIANGLE:
        tecplot_type = "FETRIANGLE"
    elif ct == vtk.VTK_QUAD:
        tecplot_type = "FEQUADRILATERAL"
        is3d = False
    elif ct == vtk.VTK_HEXAHEDRON:
        tecplot_type = "FEBRICK"
    else:
        print "WARNING: Unknown cell type in VTK unstructured grid, cannot convert"

    array_names = []
    for i in range(0,ug.GetCellData().GetNumberOfArrays()):
      array = ug.GetCellData().GetArray(i)
      array_names.append(array.GetName())

    fh.write("TITLE = " + '"' + title + '"' + "\n")
    if is3d:
        fh.write('VARIABLES = "X", "Y", "Z"')
    else:
        fh.write('VARIABLES = "X", "Y"')

    for name in array_names:
        fh.write(', "' + name + '"')
    fh.write("\n")

    points = ug.GetPoints()
    fh.write("ZONE NODES = " + str(points.GetNumberOfPoints()) + 
             ", ELEMENTS = " + str(ug.GetNumberOfCells()) + 
             ", DATAPACKING = BLOCK")
    
    if array_names:
       fh.write(", VARLOCATION = ([")
       if is3d:
           fh.write(str(4) + "-" + str(3+len(array_names)))
       else:
           fh.write(str(3) + "-" + str(2+len(array_names)))
       fh.write("] = CELLCENTERED)")

    fh.write(", ZONETYPE = " + tecplot_type + ", SOLUTIONTIME = " + str(time) + "\n")

    write_points(points, fh, 0)
    write_points(points, fh, 1)
    if is3d:
        write_points(points, fh, 2)

    for i in range(0,ug.GetCellData().GetNumberOfArrays()):
        array = ug.GetCellData().GetArray(i)
        write_array(array, fh)

    for i in range(0,ug.GetNumberOfCells()):
        pids = ug.GetCell(i).GetPointIds()        
        for j in range(0,pids.GetNumberOfIds()):
            fh.write(str(pids.GetId(j)+1) + " ")
        fh.write("\n")
    
