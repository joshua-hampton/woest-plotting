from pykml import parser


def read_kml(kml_file):
    """
    Read kml file containing WesCon grid, return co-ordinates of lines

    Args:
        kml_file (str): name of kml file with grid data

    Returns:
        tuple containing

        - outer_lines (list): co-ordinate pairs for lines describing the outer grid
        - horizontal_lines (list): co-ordinate pairs for lines describing the
              horizontal lines of the grid
        - vertical_lines (list): co-ordinate pairs for lines describing the
              vertical lines of the grid
    """
    with open(kml_file) as f:
        doc = parser.parse(f)
    outer_lines = []
    h_lines = []
    v_lines = []
    storm_boxes = []
    labels = {}
    for i in doc.iter():
        if "Document" in i.__dir__():
            for f in i.Document.Folder:
                if f.name == "GridOutline":
                    for p in f.Placemark:
                        outer_lines.append(
                            (
                                [
                                    i.strip()
                                    for i in str(p.LineString.coordinates)
                                    .strip()
                                    .split("\n")
                                ]
                            )
                        )
                elif f.name == "Verticals":
                    for p in f.Placemark:
                        v_lines.append(
                            (
                                [
                                    i.strip()
                                    for i in str(p.LineString.coordinates)
                                    .strip()
                                    .split("\n")
                                ]
                            )
                        )
                elif f.name == "Horizontals":
                    for p in f.Placemark:
                        h_lines.append(
                            (
                                [
                                    i.strip()
                                    for i in str(p.LineString.coordinates)
                                    .strip()
                                    .split("\n")
                                ]
                            )
                        )
                elif f.name == "StormBoxes":
                    for p in f.Placemark:
                        storm_boxes.append(
                            (
                                [
                                    i.strip()
                                    for i in str(p.LineString.coordinates)
                                    .strip()
                                    .split(" ")
                                ]
                            )
                        )
                elif f.name == "Labels":
                    for p in f.Placemark:
                        label_name = p.name
                        label_coord = p.Point.coordinates
                        labels[label_name] = label_coord
    return outer_lines, h_lines, v_lines, storm_boxes, labels
