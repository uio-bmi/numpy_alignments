import datetime
from numpy.random import randint

def make_report(directory, ids, names, colors):
    out = """
    <html>
    <head>
        <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.2.1/css/bootstrap.min.css" integrity="sha384-GJzZqFGwb1QTTN6wy59ffF1BuGJpLSa9DkKMp0DgiMDm4iYMj70gZWKYbI706tWS" crossorigin="anonymous">

        <style>
            .image_box{
                display: table-cell;
                padding: 0px;

        </style>
    </head>

    <body>

    <div style='width: 1800px;'>
    """

    images = [("All reads", "all"), ("Reads with variants", "variants"), ("Reads without variants", "nonvariants")]

    header = "Report %s" % str(datetime.datetime.now())
    out += "<div align='center' style='width: 1800px; font-size: 2em; font-weight: bold;'>" + header + "</div>"

    directory = "./"
    for title, image_postfix in images:
        image_url = image_postfix + ".html?r=" + str(randint(0,9999999999))
        width = 600
        if "novel" in image_postfix:
            width = 600

        out += """<div class='image_box'>
        <div align='center'><p style='font-size: 1.7em; position: relative; top: 100px;'><b>%s</b></p></div>
        <iframe style='width: %dpx; height: %dpx; border: none;' src='%s'></iframe>
        </div>
        """ % (title, width, width, image_url)


    out += "<div style='font-size: 1.5em; text-align: center; margin-top: -30px;'>"
    for id in ids:
        color = colors[id]
        out += "<span style='margin-right: 30px'><font color='" + color + "'>&#9644; " + names[id] + "</font></span>"

    out += "</div>"

    # Make a separate part for all reads
    out += """
    <div style='margin-top: 100px; margin-left: 50px;'>
    <h2>All reads in separate plot</h2>
    <iframe style='width: 600px; height: 600px; border: none; float: left;' src='all.html?r=""" + str(randint(0, 9999999999)) + """'></iframe>
    <p style='margin-top: 100px;'>
    """

    for id in ids:
        color = colors[id]
        out += "<span style='margin-right: 30px; font-size: 1.5em;'><font color='" + color + "'>&#9644; " + names[id] + "</font></span><br>"

    out += """
    </div>
    """

    out += """</div></body>
    </html>"""

    return out
