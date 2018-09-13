def field_to_csv(odb_handle, elset_handle, field_name, step_number=-1, frame_number=-1):
    """
    Outputs the given field at given frame number in odb file.

    TODO:
    - Make similar for nodal quantities
    - Test
    """
    from abaqusConstants import CENTROID
    import csv

    frame_handle = odb_handle.steps.items()[step_number][1].frames[frame_number]
    field = frame_handle.fieldOutputs[field_name].getSubset(region=elset_handle, position=CENTROID)
    file_name = 'CENTROID-{}-{}-{}-{}.csv'.format(elset_handle.name, field.name, step_number, frame_number)
    print("Exporting {}".format(file_name))
    with open(file_name, 'w') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=',')
        csv_writer.writerow(('Step time:{}'.format(frame_handle.frameValue),))
        csv_writer.writerow(('Element_id',) + field.componentLabels)
        for point in field.values:
            csv_writer.writerow((point.elementLabel,) + tuple(point.data))


def field_history_to_csv(odb_handle, elset_handle, field_name, steps=[-1], frequency=1):
    """
    Creates one csv file per frame distributed with given frequency in the given step.

    TODO:
    - Test
    - Perhaps create a folder?
    """

    for step in steps:
        step_handle = odb_handle.steps.items()[step][1]
        frames = distributed_selection(range(len(step_handle.frames)), frequency)
        #frames = [step.frames[-1]]
        for frame in frames:
            field_to_csv(odb_handle, elset_handle, field_name, step, frame)


def distributed_selection(lst, num):
    """
    Selects evenly elements from the given list with the given number of elements.
    Ensures that the last point are included. The first point is not necessary due
    to being part of the previous time step.

    For two points, the mid and last point should be selected...
    
    TODO:
    - Test
    """

    if num >= len(lst):
        return lst

    selection = []

    if num > 1:
        selection = [lst[1]] # Yes, not the first element.
        if num < len(lst):
            #q, _ = divmod(len(lst), num-1)
            q = len(lst)//(num-1)
            selection += [lst[i] for i in range(q, len(lst)-1, q)]

    selection.append(lst[-1])     

    return selection


def csv_to_legacy_vtk(csv_file_name, vtk_file_name):
    """
    Translates a csv field to vtk legacy cell data..

    TODO:
    - Support for partial fields
    """

    import csv

    with open(csv_file_name, 'r') as csv_file:
        reader = csv.reader(csv_file)
        num_rows = sum([1 for row in reader]) - 1

        if num_rows < 0:
            raise Exception("File empty!")
        
        csv_file.seek(0)
        header = next(reader)
        field_name = header[2]
        print("Header: {}".format(header))
        
        print("Number of rows: {}".format(num_rows))
    	
        with open(vtk_file_name, 'w+') as vtk_file:
            vtk_file.write('\n\nCELL_DATA {}'.format(num_rows))
            vtk_file.write('\nSCALARS {}'.format(field_name))
            vtk_file.write('\nLOOKUP_TABLE default\n')
        
            for row in reader:
                vtk_file.write('{}\n'.format(row[1]))
    return


def element_stress_history_to_csv(odb, region, output_file='', steps='', frame=-1):
    """
    Extracts centroid stress history for an element set and stores in comma separated file.

    TODO:
    - Read field once for each time step? Might be more efficient.
    """

    try:
        from abaqusConstants import CENTROID
        element_set = odb.rootAssembly.elementSets[region.upper()]
    except ImportError as ex:
        raise ImportError('Failed to import abaqus modules!')
    except KeyError as ex:
        raise IOError('Could not find element set in odb-file: ' + str(ex) + '. Element set(s) found: ' + str([el[0] for el in odb.rootAssembly.elementSets.items()]))

    if output_file == '':
        csv_path = region + '.csv'
    else:
        csv_path = output_file

    if steps == '':
        step_keys = odb.steps.keys()
    else:
        step_keys = [odb.steps.keys()[i] for i in steps]

    # For debugging
    print('Found {} elements in region {}'.format(len(element_set.elements[0]), region))
    print('Number of history points: {}'.format(len(step_keys)))
    print('Generating CSV file '+ csv_path)

    if element_set.elements[0][0].type[1] == '3':
        dimension = '3D'
    else:
        dimension = '2D'

    with open(csv_path, 'w') as csv_file:
        csv_file.write('Centroid {} stress tensor history\n'.format(dimension))

        for element in element_set.elements[0]:
            row = str(element.label) + ','
            for step in step_keys:
                field = odb.steps[step].frames[frame].fieldOutputs['S'].getSubset(region=element, position=CENTROID) # Extrapolated and averaged?
                for component in field.values[0].data:
                    row = row + str(component) + ','
            row = row + '\n'
            csv_file.write(row)

    return


