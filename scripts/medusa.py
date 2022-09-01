#!/usr/bin/env python3


import os
import sys
import math
import datetime
import argparse
import h5py
import pandas as pd
import numpy as np
import functools

from tensorflow.keras.models import model_from_json
from tensorflow.keras.utils import to_categorical
from tensorflow.keras.losses import CategoricalCrossentropy
from tensorflow.keras.optimizers import Adam

from itertools import product
from sklearn.metrics import matthews_corrcoef


# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("-i", help="Path to the merge feature vector file.", type=str, required=True)
parser.add_argument("-o", "--output", help="Path to the output directory.", type=str, required=True)
parser.add_argument("-m", "--models", help="Path to the models directory.", type=str, required=True)
parser.add_argument("-f", "--fasta", help="Path to the fasta file.", type=str, required=True)
parser.add_argument("-S", "--strict", help="Strict prediction.", action="store_true", default=True)
parser.add_argument("-NS", "--nonstrict", help="Non strict prediction.", action="store_true", default=True)
parser.add_argument("-3", "--three", help="3 classes prediction.", action="store_true", default=True)
parser.add_argument("-5", "--five", help="5 classes prediction.", action="store_true", default=True)
parser.add_argument("-d", "--dim", metavar="D", type=int, nargs=2, help="Rows Columns", required=True)
args = parser.parse_args()

MODEL = args.models

#Check if the merge vector file is valid
merge = args.i
if not os.path.isfile(merge):
    sys.exit(f"{merge} does not exist.\n"
			  "Please enter a valid merge vector file.")

OUTPUT = args.output

fasta = args.fasta
if not os.path.isfile(fasta):
    sys.exit(f"{fasta} does not exist.\n"
			  "Please enter a valid fasta file.")

if not (args.strict or args.nonstrict or three or five):
    parser.error("No model selected.")

model_list = []
if args.strict:
    model_list.append("S")
if args.nonstrict:
    model_list.append("NS")
if args.three:
    model_list.append("3")
if args.five:
    model_list.append("5")


ROWS, COLS = args.dim


if __name__ == "__main__":

    #Get entry files path
    ENTRY = merge

    ########################### LOADING MODELS #################################
    if "S" in model_list:

        PATH_MODEL_JSON = os.path.join(MODEL, "S", "model.json")
        PATH_MODEL_H5 = os.path.join(MODEL, "S", "weights_S.h5")

        # Load JSON arch
        with open(PATH_MODEL_JSON, "r") as json_file:
            loaded_model_json = json_file.read()
        model_S = model_from_json(loaded_model_json)

        # Load weights from H5
        model_S.load_weights(PATH_MODEL_H5)

        # Compile model
        model_S.compile(loss=CategoricalCrossentropy(),
                        optimizer=Adam())

    if "NS" in model_list:

        PATH_MODEL_JSON = os.path.join(MODEL, "NS", "model.json")
        PATH_MODEL_H5 = os.path.join(MODEL, "NS", "weights_NS.h5")

        # Load JSON arch
        with open(PATH_MODEL_JSON, "r") as json_file:
            loaded_model_json = json_file.read()
        model_NS = model_from_json(loaded_model_json)

        # Load weights from H5
        model_NS.load_weights(PATH_MODEL_H5)

        # Compile model
        model_NS.compile(loss=CategoricalCrossentropy(),
                        optimizer=Adam())

    if "3" in model_list:

        PATH_MODEL_JSON = os.path.join(MODEL, "3", "model.json")
        PATH_MODEL_H5 = os.path.join(MODEL, "3", "weights_3.h5")

        # Load JSON arch
        with open(PATH_MODEL_JSON, "r") as json_file:
            loaded_model_json = json_file.read()
        model_3 = model_from_json(loaded_model_json)

        # Load weights from H5
        model_3.load_weights(PATH_MODEL_H5)

        # Compile model
        model_3.compile(loss=CategoricalCrossentropy(),
                        optimizer=Adam())

    if "5" in model_list:

        PATH_MODEL_JSON = os.path.join(MODEL, "5", "model.json")
        PATH_MODEL_H5 = os.path.join(MODEL, "5", "weights_5.h5")

        # Load JSON arch
        with open(PATH_MODEL_JSON, "r") as json_file:
            loaded_model_json = json_file.read()
        model_5 = model_from_json(loaded_model_json)

        # Load weights from H5
        model_5.load_weights(PATH_MODEL_H5)

        # Compile model
        model_5.compile(loss=CategoricalCrossentropy(),
                        optimizer=Adam())




    #Load X
    X = np.loadtxt(ENTRY)
    x_test = np.reshape(X, (len(X), ROWS, COLS))

    #Load fasta file
    fasta_file = np.loadtxt(fasta, dtype="str", skiprows=1)
    seq_AA = fasta_file.tolist()

    #Create output directory
    os.makedirs(OUTPUT, exist_ok=False)


    if "S" in model_list:
        y_pred = model_S.predict(x_test, verbose=0)
        y_pred_classes = np.argmax(y_pred, axis=1)
        prob_y_pred = np.amax(y_pred, axis=1)

        # Create a dataframe from zipped list
        zippedList =  list(zip(seq_AA, y_pred_classes, prob_y_pred, y_pred[:, 0], y_pred[:, 1]))
        df = pd.DataFrame(zippedList, columns = ["res", "S", "P_max", "P_0", "P_1"])
        df.to_csv(os.path.join(OUTPUT, "S_prediction.csv"), index=False, float_format="%.2f", sep="\t")


    if "NS" in model_list:
        y_pred = model_NS.predict(x_test, verbose=0)
        y_pred_classes = np.argmax(y_pred, axis=1)
        prob_y_pred = np.amax(y_pred, axis=1)

        # Create a dataframe from zipped list
        zippedList =  list(zip(seq_AA, y_pred_classes, prob_y_pred, y_pred[:, 0], y_pred[:, 1]))
        df = pd.DataFrame(zippedList, columns = ["res", "NS", "P_max", "P_0", "P_1"])
        df.to_csv(os.path.join(OUTPUT, "NS_prediction.csv"), index=False, float_format="%.2f", sep="\t")


    if "3" in model_list:
        y_pred = model_3.predict(x_test, verbose=0)
        y_pred_classes = np.argmax(y_pred, axis=1)
        prob_y_pred = np.amax(y_pred, axis=1)

        # Create a dataframe from zipped list
        zippedList =  list(zip(seq_AA, y_pred_classes, prob_y_pred, y_pred[:, 0], y_pred[:, 1], y_pred[:, 2]))
        df = pd.DataFrame(zippedList, columns = ["res", "pred_3", "P_max", "P_0", "P_1", "P_2"])
        df.to_csv(os.path.join(OUTPUT, "3_prediction.csv"), index=False, float_format="%.2f", sep="\t")


    if "5" in model_list:
        y_pred = model_5.predict(x_test, verbose=0)
        y_pred_classes = np.argmax(y_pred, axis=1)
        prob_y_pred = np.amax(y_pred, axis=1)

        # Create a dataframe from zipped list
        zippedList =  list(zip(seq_AA, y_pred_classes, prob_y_pred, y_pred[:, 0], y_pred[:, 1], y_pred[:, 2], y_pred[:, 3], y_pred[:, 4]))
        df = pd.DataFrame(zippedList, columns = ["res", "pred_5", "P_max", "P_0", "P_1", "P_2", "P_3", "P_4"])
        df.to_csv(os.path.join(OUTPUT, "5_prediction.csv"), index=False, float_format="%.2f", sep="\t")

header = ""
with open(fasta, "r") as f:
    header = f.readline().strip()[1:]

os.makedirs(f"{OUTPUT}/../html", exist_ok=False)

with open(f"{OUTPUT}/../html/results.html", "w") as f:
    f.write(f"""
	<!DOCTYPE html>
	<html lang="en">

	<head>
        <meta charset="utf-8" />
        <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no" />
        <meta name="description" content="" />
        <meta name="author" content="" />

        <title>MEDUSA | Results</title>

        <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.5.0/css/bootstrap.min.css" integrity="sha384-9aIt2nRpC12Uk9gS9baDl411NQApFmC26EwAOH8WgZl5MYYxFfc+NcPb1dKGj7Sk" crossorigin="anonymous"/>
        <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.4.1/jquery.min.js"></script>
        <script src="https://cdn.jsdelivr.net/npm/popper.js@1.16.0/dist/umd/popper.min.js" integrity="sha384-Q6E9RHvbIyZFJoft+2mJbHaEWldlvI9IOYy5n3zV9zzTtmI3UksdQRVvoxMfooAo" crossorigin="anonymous"></script>
        <script src="https://stackpath.bootstrapcdn.com/bootstrap/4.5.0/js/bootstrap.min.js" integrity="sha384-OgVRvuATP1z7JjHLkuOU7Xw704+h835Lr+6QL9UvYjZE3Ipu6Tp75j7Bh/kR0JKI" crossorigin="anonymous"></script>
        <script defer src="https://use.fontawesome.com/releases/v5.13.0/js/all.js" integrity="sha384-ujbKXb9V3HdK7jcWL6kHL1c+2Lj4MR4Gkjl7UtwpSHg/ClpViddK9TI7yU53frPN" crossorigin="anonymous"></script>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/Chart.js/2.9.3/Chart.min.js"></script>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/d3/5.16.0/d3.min.js"></script>
        <script src="https://cdn.jsdelivr.net/npm/html2canvas@0.5.0-beta4/dist/html2canvas.min.js" integrity="sha256-jP/zNgDdOLDvHgG/IPyrDDazrf7b+P6CMMZ/S3EZ6zc=" crossorigin="anonymous"></script>

        <!-- Bootstrap core CSS -->
        <link href="https://fonts.googleapis.com/css2?family=Overpass:ital,wght@0,100;0,200;0,300;0,400;0,600;0,700;0,800;0,900;1,100;1,200;1,300;1,400;1,600;1,700;1,800;1,900&display=swap" rel="stylesheet"/>
        <link href="https://fonts.googleapis.com/css2?family=Overpass:ital,wght@0,100;0,200;0,300;0,400;0,600;0,700;0,800;0,900;1,100;1,200;1,300;1,400;1,600;1,700;1,800;1,900&family=Roboto:ital,wght@0,100;0,300;0,400;0,500;0,700;0,900;1,100;1,300;1,400;1,500;1,700;1,900&display=swap" rel="stylesheet"/>
        <link href="https://fonts.googleapis.com/css2?family=Roboto+Mono:ital,wght@0,100;0,200;0,300;0,400;0,500;0,600;0,700;1,100;1,200;1,300;1,400;1,500;1,600;1,700&display=swap" rel="stylesheet"/>

        <!-- Custom styles for this template -->
        <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css" />
        <link href="../../../html/css/custom_features.css" rel="stylesheet" />
    </head>


        <nav class="navbar navbar-expand-lg navbar-light bg-white fixed-top mb-4 pb-0">
            <div class="container shadow-sm" id="nav-container">
                <a class="navbar-brand"><img src="../../../html/images/MEDUSA_logo.png" alt=""></a>
            </div>
        </nav>
	<body>

		<!-- Page Content -->
		<div class="container mt-10" id="principal_container">
			<div class="alert alert-success" role="alert">
			    <h1 class="alert-heading">Results</h1>
                            <p class="text-break m-0"><strong>Query name: </strong> {header}</p>
                            <p class="text-break m-0"><strong>Query sequence: </strong> {"".join(seq_AA)}</p>
	                </div>
		<nav id="res-nav">
            <div class="nav nav-tabs nav-fill mt-6 mb-5" id="nav-tab" role="tablist">
                <a class="nav-item nav-link active" id="nav-s-tab" data-toggle="tab" href="#nav-tabS" role="tab" aria-controls="nav-tabS" aria-selected="true">Strict threshold</a>
                <a class="nav-item nav-link" id="nav-ns-tab" data-toggle="tab" href="#nav-tabNS" role="tab" aria-controls="nav-tabNS" aria-selected="false">Non-strict threshold</a>
                <a class="nav-item nav-link" id="nav-3pred-tab" data-toggle="tab" href="#nav-tab3pred" role="tab" aria-controls="nav-tab3pred" aria-selected="false">3 classes prediction</a>
                <a class="nav-item nav-link" id="nav-5pred-tab" data-toggle="tab" href="#nav-tab5pred" role="tab" aria-controls="nav-tab5pred" aria-selected="false">5 classes prediction</a>
            </div>
        </nav>
		<div class="tab-content">
			<div class="tab-pane fade show active" id="nav-tabS" role="tabpanel" aria-labelledby="nav-tabS">
				<h1 class="alert-heading">Sequence plot</h1>
                <a href="javascript:void(0)" class="btn btn-outline-primary btn-sm float-right mt-4" id="btn-download-chartS">Download PNG</a>
                <div class="d-flex justify-content-start">
                    <p class="d-flex align-items-center">Confidence of prediction: </p>
                    <img src="../../../html/images/legend-2cl.png" alt="legend-2cl" width="150" height="60">
                </div>
                <div id="divChartS">
                </div>
                <h1 class="mt-6 alert-heading">Flexibility class distribution</h1>
				<a href="javascript:void(0)" class="btn btn-outline-primary btn-sm float-right mt-4" id="btn-download-pieChartS">Download PNG</a>
                <div id="divPieChartS" style="width: 60%">
                    <canvas id="pie_chart_S" height="150"></canvas>
                    <div  class="d-flex justify-content-left">
                        <img src="../../../html/images/legend_2cl_piechart.png" alt="legend_2cl_piechart" style="height: 70px" class="ml-5"/>
                    </div>
                </div>
			</div>
			<div class="tab-pane fade" id="nav-tabNS" role="tabpanel" aria-labelledby="nav-tabNS">
				<h1 class="alert-heading">Sequence plot</h1>
                <a href="javascript:void(0)" class="btn btn-outline-primary btn-sm float-right mt-4" id="btn-download-chartNS">Download PNG</a>
                <div class="d-flex justify-content-start">
                    <p class="d-flex align-items-center">Confidence of prediction: </p>
                    <img src="../../../html/images/legend-2cl.png" alt="legend-2cl" width = "150" height="60">
                </div>
                <div id="divChartNS">
                </div>
                <h1 class="mt-6 alert-heading">Flexibility class distribution</h1>
				<a href="javascript:void(0)" class="btn btn-outline-primary btn-sm float-right mt-4" id="btn-download-pieChartNS">Download PNG</a>
                <div id="divPieChartNS" style="width: 60%">
                    <canvas id="pie_chart_NS" height="150"></canvas>
                    <div  class="d-flex justify-content-left">
                        <img src="../../../html/images/legend_2cl_piechart.png" alt="legend_2cl_piechart" style="height: 70px" class="ml-5"/>
                    </div>
                </div>
			</div>
			<div class="tab-pane fade" id="nav-tab3pred" role="tabpanel" aria-labelledby="nav-tab3pred">
				<h1 class="alert-heading">Sequence plot</h1>
                <a href="javascript:void(0)" class="btn btn-outline-primary btn-sm float-right mt-4" id="btn-download-chart3pred">Download PNG</a>
                <div class="d-flex justify-content-start">
                    <p class="d-flex align-items-center">Confidence of prediction: </p>
                    <img src="../../../html/images/legend-3cl.png" alt="legend-3cl" width = "150" height="60">
                </div>
                <div id="divChart3pred">
                </div>
                <h1 class="mt-6 alert-heading">Flexibility class distribution</h1>
				<a href="javascript:void(0)" class="btn btn-outline-primary btn-sm float-right mt-4" id="btn-download-pieChart3pred">Download PNG</a>
                <div id="divPieChart3pred" style="width: 60%">
                    <canvas id="pie_chart_3pred" height="150"></canvas>
                    <div  class="d-flex justify-content-left">
                        <img src="../../../html/images/legend_3cl_piechart.png" alt="legend_3cl_piechart" style="height: 70px" class="ml-5"/>
                    </div>
                </div>
			</div>
			<div class="tab-pane fade" id="nav-tab5pred" role="tabpanel" aria-labelledby="nav-tab5pred">
				<h1 class="alert-heading">Sequence plot</h1>
                <a href="javascript:void(0)" class="btn btn-outline-primary btn-sm float-right mt-4" id="btn-download-chart5pred">Download PNG</a>
                <div class="d-flex justify-content-start">
                    <p class="d-flex align-items-center">Confidence of prediction: </p>
                    <img src="../../../html/images/legend-5cl.png" alt="legend-5cl" width = "150" height="60">
                </div>
                <div id="divChart5pred">
                </div>
                <h1 class="mt-6 alert-heading">Flexibility class distribution</h1>
				<a href="javascript:void(0)" class="btn btn-outline-primary btn-sm float-right mt-4" id="btn-download-pieChart5pred">Download PNG</a>
                <div id="divPieChart5pred" style="width: 60%">
                    <canvas id="pie_chart_5pred" height="150"></canvas>
                    <div  class="d-flex justify-content-left">
                        <img src="../../../html/images/flexible_rigid_legend_arrow.png" alt="legend_5cl_piechart" style="height: 70px" class="ml-5"/>
                    </div>
                </div>
			</div>
		</div>
		</div>

		<script>

			// When page is scrolled towards bottom, a line appears
			// to create separation with nav bar
			//
			$(window).scroll(function() {{
				var scroll = $(window).scrollTop();

				 //>=, not <=
				if (scroll >= 100) {{
					$("#nav-container").addClass("shadow-sm");
				}} else{{
					$("#nav-container").removeClass("shadow-sm");
				}}
			}});

			// Create all "Download PNG" buttons of sequence plot chart
			//
			$("#btn-download-chartS").click(function() {{
                html2canvas($("#divChartS"), {{
                  onrendered: function(canvas) {{
                    saveAs(canvas.toDataURL(), 'strict_threshold_pred_chart.png');
                  }}
                }});
              }});

            $("#btn-download-chartNS").click(function() {{
                html2canvas($("#divChartNS"), {{
                  onrendered: function(canvas) {{
                    saveAs(canvas.toDataURL(), 'nonstrict_threshold_pred_chart.png');
                  }}
                }});
              }});

            $("#btn-download-chart3pred").click(function() {{
                html2canvas($("#divChart3pred"), {{
                  onrendered: function(canvas) {{
                    saveAs(canvas.toDataURL(), '3_classes_pred_chart.png');
                  }}
                }});
              }});

            $("#btn-download-chart5pred").click(function() {{
                html2canvas($("#divChart5pred"), {{
                  onrendered: function(canvas) {{
                    saveAs(canvas.toDataURL(), '5_classes_pred_chart.png');
                  }}
                }});
              }});


			// Create all "Download PNG" buttons of Flexibility class distribution pie chart
            //
			$("#btn-download-pieChartS").click(function() {{
                html2canvas($("#divPieChartS"), {{
                  onrendered: function(canvas) {{
                    saveAs(canvas.toDataURL(), 'strict_threshold_pie_chart.png');
                  }}
                }});
              }});

            $("#btn-download-pieChartNS").click(function() {{
                html2canvas($("#divPieChartNS"), {{
                  onrendered: function(canvas) {{
                    saveAs(canvas.toDataURL(), 'nonstrict_threshold_pie_chart.png');
                  }}
                }});
              }});

            $("#btn-download-pieChart3pred").click(function() {{
                html2canvas($("#divPieChart3pred"), {{
                  onrendered: function(canvas) {{
                    saveAs(canvas.toDataURL(), '3_classes_pie_chart.png');
                  }}
                }});
              }});

			$("#btn-download-pieChart5pred").click(function() {{
                html2canvas($("#divPieChart5pred"), {{
                  onrendered: function(canvas) {{
                    saveAs(canvas.toDataURL(), '5_classes_pie_chart.png');
                  }}
                }});
              }});




			// Create pseudo <a> link to download PNG image
			function saveAs(uri, filename) {{
				var link = document.createElement('a');
				if (typeof link.download === 'string') {{
					link.href = uri;
					link.download = filename;

					//Firefox requires the link to be in the body
					document.body.appendChild(link);

					//simulate click
					link.click();

					//remove the link when done
					document.body.removeChild(link);
				}} else {{
					window.open(uri);
				}}
			}}






			// PIE CHART
d3.tsv("../prediction/S_prediction.csv").then(

				function(csv_data) {{
					var flex = csv_data.map(function(d) {{ return d.S; }});
					// parse class column, cast string into numerical value
					// and calculate occurences of each classes for pie chart
					var values = {{}};
					for (var i=0; i < flex.length; i++) {{
						values[parseInt(flex[i], 10)] = (values[parseInt(flex[i], 10)] || 0) + 1;
					}};
					var theHelp = Chart.helpers;
					var chart_S = new Chart("pie_chart_S", {{
						type: "pie",
						data: {{
						  labels: ["Rigid", "Flexible"],
						  datasets: [
										{{
											data: Object.values(values),
											backgroundColor: ["rgba(0, 78, 137, 1)", "rgba(234, 115, 23, 1)"],
										}}
									]
						}},
						options: {{
                            legend: {{
                                display: true,
                                position: "right",
                                labels: {{
                                    generateLabels: function(chart) {{
                                        var data = chart.data;
                                        if (data.labels.length && data.datasets.length) {{
                                          return data.labels.map(function(label, i) {{
                                            var meta = chart.getDatasetMeta(0);
                                            var ds = data.datasets[0];
                                            var sum = ds["data"].reduce((pv, cv) => pv + cv, 0);
                                            var arc = meta.data[i];
                                            var custom = arc && arc.custom || {{}};
                                            var getValueAtIndexOrDefault = theHelp.getValueAtIndexOrDefault;
                                            var arcOpts = chart.options.elements.arc;
                                            var fill = custom.backgroundColor ? custom.backgroundColor : getValueAtIndexOrDefault(ds.backgroundColor, i, arcOpts.backgroundColor);
                                            var stroke = custom.borderColor ? custom.borderColor : getValueAtIndexOrDefault(ds.borderColor, i, arcOpts.borderColor);
                                            var bw = custom.borderWidth ? custom.borderWidth : getValueAtIndexOrDefault(ds.borderWidth, i, arcOpts.borderWidth);
                                              return {{
                                              // And finally :
                                              text: ((ds.data[i]*100)/sum).toFixed(0) + "% of class " + label,
                                              fillStyle: fill,
                                              strokeStyle: stroke,
                                              lineWidth: bw,
                                              hidden: isNaN(ds.data[i]) || meta.data[i].hidden,
                                              index: i
                                            }};
                                          }});
                                        }}
                                        return [];
                                    }}
                                }}
                            }},
                            title: {{
                                display: true,
                                fontFamily: "Roboto",
                                fontSize: 16,
                                fontColor: "black",
                                fontStyle: "normal",
                                lineHeight: 1.5,
                                text: ["Global flexibility of the target sequence:", "residues are grouped into 2 predicted classes of flexibility"]
                            }},
							tooltips: {{
							callbacks: {{
								label: function(tooltipItem, data) {{
									return " "+values[tooltipItem.index]+" residues";
								}}
							  }}

						}}
						}}
					}});
				}}
			);

			// POINT MAP FOR S_PREDICTIONS
d3.tsv("../prediction/S_prediction.csv").then(

				function(csv_data) {{
					var residues = csv_data.map(function(d) {{ return d.res; }});
					var flex = csv_data.map(function(d) {{ return d.S; }});
					var pred = csv_data.map(function(d) {{ return d.P_max; }});
					var numberOfCharts = Math.ceil(residues.length / 50);
                    var idx = [];
                    for (var i=1; i <= residues.length; i++) {{
                        idx.push(i);
                    }}
					var yLabels = {{
							0: "Rigid",
							1: "Flexible"
						}};

					for (var i=1; i<=numberOfCharts; i++) {{
						var canvas = document.createElement("canvas");
						canvas.id = "canvasChartS_"+i;
						canvas.height = 90;
						canvas.style = "margin-top: 25px";
						canvas.width = document.getElementById("principal_container").offsetWidth - 30;
						var mutate1_2 = false,
                            mutate3 = false,
                            mutate4 = false;
						if ((residues.length > 50) && (i == numberOfCharts) && ((residues.length % 50)!=0)) {{
								if (residues.length % 50 >= 5) {{
									canvas.width = (document.getElementById("principal_container").offsetWidth - 30)*(residues.length % 50)*1.0 / 50 + 75*(1- (residues.length % 50)*1.0/50) ;
								}} else if (residues.length % 50 == 4) {{
                                    canvas.width = 135;
                                    canvas.style = "margin-left: 0.3rem !important";
                                    var mutate4 = true;
                                }} else if (residues.length % 50 == 3) {{
                                    canvas.width = 135;
                                    canvas.style = "margin-left: -0.7rem !important";
                                    var mutate3 = true;
                                }} else if (residues.length % 50 == 2 || residues.length % 50 == 1) {{
                                    canvas.width = 135;
                                    canvas.style = "margin-left: -1.7rem !important";
                                    var mutate1_2 = true;
                                }}
						}}
						div = document.getElementById("divChartS");
						div.appendChild(canvas);

						var chart = new Chart("canvasChartS_"+i, {{
							type: "line",
							data: {{
							  labels: residues.slice((i*50)-50, (i*50)),
							  datasets: [
											{{
												data: flex.slice((i*50)-50, (i*50)),
												flex: pred.slice((i*50)-50, (i*50)),
                                                idx: idx.slice((i*50)-50, (i*50)),
												tension: 0,
												fill: false,
												borderWidth: 0,
												pointBorderWidth: 0,
												pointRadius: 10,
												pointHoverRadius: 11,
												showLine: false,
												pointBackgroundColor: []
											}}
										]
							}},
							options: {{
								legend: {{
									display: false
								}},
								responsive: false,
								elements: {{
									point: {{
										pointStyle: "rect",
										borderWidth: 0,
										hoverBorderWidth: 0
									}}
								}},
								scales: {{
											xAxes: [{{
												ticks: {{
													autoSkip: false,
													maxRotation: 0,
													fontStyle: "bold",
													fontSize: 15,
													callback: function(value, index, values) {{
														if ((index+1) % 10 == 0 || (index+1) % 50 == 1) {{
															return [value,index+((i-1)*50)+1];
														}}
														else {{
															return value;
														}}
													}}
												}}
											}}],
											yAxes: [{{
												ticks: {{
													beginAtZero: true,
													fontSize: 16,
													padding: 10,
													callback: function(value, index, values) {{
														return yLabels[value];
													}}
												}}
											}}]
										}},
								tooltips: {{
									callbacks: {{
										title: function(tooltipItems, data) {{
                                            return data.labels[tooltipItems[0].index] + " (" + data.datasets[0].idx[tooltipItems[0].index] + ")";
                                        }},
                                        label: function(tooltipItem, data) {{
                                            return "Probabillity: "+((data.datasets[0].flex[tooltipItem.index]*100).toFixed(0))+" %";
                                        }}
									  }}

								}}
							}}
						}});

						if (mutate4){{
                            chart.options.layout.padding.left = 21;
                            chart.update();
                        }} else if (mutate3) {{
                            chart.options.layout.padding.left = 53;
                            chart.update();
                        }} else if (mutate1_2) {{
                            chart.options.layout.padding.left = 87;
                            chart.update();
                        }}

						var chartColors = {{
							color1: "rgba(35, 98, 118, 1)",
							color2: "rgba(0, 153, 159, 1)",
							color3: "rgba(0, 209, 166, 1)",
						}};


						var dataset = chart.data.datasets[0];
						for (var j = 0; j < dataset.data.length; j++) {{
							k = j + ((i-1)*50);
							if (pred[k] < 0.6) {{
								dataset.pointBackgroundColor[j] = chartColors.color3;
							}}
							else if ((pred[k] >= 0.6) && (pred[k] <= 0.75)){{
								dataset.pointBackgroundColor[j] = chartColors.color2;
							}}
							else {{
								dataset.pointBackgroundColor[j] = chartColors.color1;
							}}
						}}
						chart.update();

						window.chart = chart;
				}};
			}});


			// PIE CHART
d3.tsv("../prediction/NS_prediction.csv").then(

				function(csv_data) {{
					var flex = csv_data.map(function(d) {{ return d.NS; }});
					// parse class column, cast string into numerical value
					// and calculate occurences of each classes for pie chart
					var values = {{}};
					for (var i=0; i < flex.length; i++) {{
						values[parseInt(flex[i], 10)] = (values[parseInt(flex[i], 10)] || 0) + 1;
					}};
					var theHelp = Chart.helpers;
					var chart_NS = new Chart("pie_chart_NS", {{
						type: "pie",
						data: {{
						  labels: ["Rigid", "Flexible"],
						  datasets: [
										{{
											data: Object.values(values),
											backgroundColor: ["rgba(0, 78, 137, 1)", "rgba(234, 115, 23, 1)"],
										}}
									]
						}},
						options: {{
                            legend: {{
                                display: true,
                                position: "right",
                                labels: {{
                                    generateLabels: function(chart) {{
                                        var data = chart.data;
                                        if (data.labels.length && data.datasets.length) {{
                                          return data.labels.map(function(label, i) {{
                                            var meta = chart.getDatasetMeta(0);
                                            var ds = data.datasets[0];
                                            var sum = ds["data"].reduce((pv, cv) => pv + cv, 0);
                                            var arc = meta.data[i];
                                            var custom = arc && arc.custom || {{}};
                                            var getValueAtIndexOrDefault = theHelp.getValueAtIndexOrDefault;
                                            var arcOpts = chart.options.elements.arc;
                                            var fill = custom.backgroundColor ? custom.backgroundColor : getValueAtIndexOrDefault(ds.backgroundColor, i, arcOpts.backgroundColor);
                                            var stroke = custom.borderColor ? custom.borderColor : getValueAtIndexOrDefault(ds.borderColor, i, arcOpts.borderColor);
                                            var bw = custom.borderWidth ? custom.borderWidth : getValueAtIndexOrDefault(ds.borderWidth, i, arcOpts.borderWidth);
                                              return {{
                                              // And finally :
                                              text: ((ds.data[i]*100)/sum).toFixed(0) + "% of class " + label,
                                              fillStyle: fill,
                                              strokeStyle: stroke,
                                              lineWidth: bw,
                                              hidden: isNaN(ds.data[i]) || meta.data[i].hidden,
                                              index: i
                                            }};
                                          }});
                                        }}
                                        return [];
                                    }}
                                }}
                            }},
                            title: {{
                                display: true,
                                fontFamily: "Roboto",
                                fontSize: 16,
                                fontColor: "black",
                                fontStyle: "normal",
                                lineHeight: 1.5,
                                text: ["Global flexibility of the target sequence:", "residues are grouped into 2 predicted classes of flexibility"]
                            }},
							tooltips: {{
							callbacks: {{
								label: function(tooltipItem, data) {{
									return " "+values[tooltipItem.index]+" residues";
								}}
							  }}

						}}
						}}
					}});
				}}
			);

			// POINT MAP FOR NS_PREDICTIONS
d3.tsv("../prediction/NS_prediction.csv").then(

				function(csv_data) {{
					var residues = csv_data.map(function(d) {{ return d.res; }});
					var flex = csv_data.map(function(d) {{ return d.NS; }});
					var pred = csv_data.map(function(d) {{ return d.P_max; }});
					var numberOfCharts = Math.ceil(residues.length / 50);
                    var idx = [];
                    for (var i=1; i <= residues.length; i++) {{
                        idx.push(i);
                    }}
					var yLabels = {{
							0: "Rigid",
							1: "Flexible"
						}};

					for (var i=1; i<=numberOfCharts; i++) {{
						var canvas = document.createElement("canvas");
						canvas.id = "canvasChartNS_"+i;
						canvas.height = 90;
						canvas.style = "margin-top: 25px";
						canvas.width = document.getElementById("principal_container").offsetWidth - 30;
						var mutate1_2 = false,
                            mutate3 = false,
                            mutate4 = false;
						if (residues.length > 50 && i == numberOfCharts && ((residues.length % 50)!=0)) {{
								if (residues.length % 50 >= 5) {{
									canvas.width = (document.getElementById("principal_container").offsetWidth - 30)*(residues.length % 50)*1.0 / 50 + 75*(1- (residues.length % 50)*1.0/50);
								}} else if (residues.length % 50 == 4) {{
                                    canvas.width = 135;
                                    canvas.style = "margin-left: 0.3rem !important";
                                    var mutate4 = true;
                                }} else if (residues.length % 50 == 3) {{
                                    canvas.width = 135;
                                    canvas.style = "margin-left: -0.7rem !important";
                                    var mutate3 = true;
                                }} else if (residues.length % 50 == 2 || residues.length % 50 == 1) {{
                                    canvas.width = 135;
                                    canvas.style = "margin-left: -1.7rem !important";
                                    var mutate1_2 = true;
                                }}
						}}

						div = document.getElementById("divChartNS");
						div.appendChild(canvas);

						var chart = new Chart("canvasChartNS_"+i, {{
							type: "line",
							data: {{
							  labels: residues.slice((i*50)-50, (i*50)),
							  datasets: [
											{{
												data: flex.slice((i*50)-50, (i*50)),
												flex: pred.slice((i*50)-50, (i*50)),
                                                idx: idx.slice((i*50)-50, (i*50)),
												tension: 0,
												fill: false,
												borderWidth: 0,
												pointBorderWidth: 0,
												pointRadius: 10,
												pointHoverRadius: 11,
												showLine: false,
												pointBackgroundColor: []
											}}
										]
							}},
							options: {{
								legend: {{
									display: false
								}},
								responsive: false,
								elements: {{
									point: {{
										pointStyle: "rect",
										borderWidth: 0,
										hoverBorderWidth: 0
									}}
								}},
								scales: {{
											xAxes: [{{
												ticks: {{
													autoSkip: false,
													maxRotation: 0,
													fontStyle: "bold",
													fontSize: 15,
													callback: function(value, index, values) {{
														if ((index+1) % 10 == 0 || (index+1) % 50 == 1) {{
															return [value,index+((i-1)*50)+1];
														}}
														else {{
															return value;
														}}
													}}
												}}
											}}],
											yAxes: [{{
												ticks: {{
													beginAtZero: true,
													fontSize: 16,
													padding: 10,
													callback: function(value, index, values) {{
														return yLabels[value];
													}}
												}}
											}}]
										}},
								tooltips: {{
									callbacks: {{
										title: function(tooltipItems, data) {{
                                            return data.labels[tooltipItems[0].index] + " (" + data.datasets[0].idx[tooltipItems[0].index] + ")";
                                        }},
                                        label: function(tooltipItem, data) {{
                                            return "Probabillity: "+((data.datasets[0].flex[tooltipItem.index]*100).toFixed(0))+" %";
                                        }}
									  }}

								}}
							}}
						}});

						if (mutate4){{
                            chart.options.layout.padding.left = 21;
                            chart.update();
                        }} else if (mutate3) {{
                            chart.options.layout.padding.left = 53;
                            chart.update();
                        }} else if (mutate1_2) {{
                            chart.options.layout.padding.left = 87;
                            chart.update();
                        }}


						var chartColors = {{
							color1: "rgba(35, 98, 118, 1)",
							color2: "rgba(0, 153, 159, 1)",
							color3: "rgba(0, 209, 166, 1)",
						}};


						var dataset = chart.data.datasets[0];
						for (var j = 0; j < dataset.data.length; j++) {{
							k = j + ((i-1)*50);
							if (pred[k] < 0.6) {{
								dataset.pointBackgroundColor[j] = chartColors.color3;
							}}
							else if ((pred[k] >= 0.6) && (pred[k] <= 0.75)){{
								dataset.pointBackgroundColor[j] = chartColors.color2;
							}}
							else {{
								dataset.pointBackgroundColor[j] = chartColors.color1;
							}}
						}}
						chart.update();

						window.chart = chart;
				}};
			}});


			// PIE CHART 3PRED
d3.tsv("../prediction/3_prediction.csv").then(

				function(csv_data) {{
					var flex = csv_data.map(function(d) {{ return d.pred_3; }});
					// parse class column, cast string into numerical value
					// and calculate occurences of each classes for pie chart
					var values = {{0: 0, 1: 0, 2: 0}};
					for (var i=0; i < flex.length; i++) {{
						values[parseInt(flex[i], 10)] = (values[parseInt(flex[i], 10)] || 0) + 1;
					}};
					var theHelp = Chart.helpers;
					var chart_3pred = new Chart("pie_chart_3pred", {{
						type: "pie",
						data: {{
						  labels: ["0", "1", "2"],
						  datasets: [
										{{
											data: Object.values(values),
											backgroundColor: ["rgba(0, 78, 137, 1)", "rgba(61, 165, 217, 1)", "rgba(234, 115, 23, 1)"],
										}}
									]
						}},
						options: {{
                            legend: {{
                                display: true,
                                position: "right",
                                labels: {{
                                    generateLabels: function(chart) {{
                                        var data = chart.data;
                                        if (data.labels.length && data.datasets.length) {{
                                          return data.labels.map(function(label, i) {{
                                            var meta = chart.getDatasetMeta(0);
                                            var ds = data.datasets[0];
                                            var sum = ds["data"].reduce((pv, cv) => pv + cv, 0);
                                            var arc = meta.data[i];
                                            var custom = arc && arc.custom || {{}};
                                            var getValueAtIndexOrDefault = theHelp.getValueAtIndexOrDefault;
                                            var arcOpts = chart.options.elements.arc;
                                            var fill = custom.backgroundColor ? custom.backgroundColor : getValueAtIndexOrDefault(ds.backgroundColor, i, arcOpts.backgroundColor);
                                            var stroke = custom.borderColor ? custom.borderColor : getValueAtIndexOrDefault(ds.borderColor, i, arcOpts.borderColor);
                                            var bw = custom.borderWidth ? custom.borderWidth : getValueAtIndexOrDefault(ds.borderWidth, i, arcOpts.borderWidth);
                                              return {{
                                              // And finally :
                                              text: ((ds.data[i]*100)/sum).toFixed(0) + "% of class " + label,
                                              fillStyle: fill,
                                              strokeStyle: stroke,
                                              lineWidth: bw,
                                              hidden: isNaN(ds.data[i]) || meta.data[i].hidden,
                                              index: i
                                            }};
                                          }});
                                        }}
                                        return [];
                                    }}
                                }}
                            }},
                            title: {{
                                display: true,
                                fontFamily: "Roboto",
                                fontSize: 16,
                                fontColor: "black",
                                fontStyle: "normal",
                                lineHeight: 1.5,
                                text: ["Global flexibility of the target sequence:", "residues are grouped into 3 predicted classes of flexibility"]
                            }},
							tooltips: {{
							callbacks: {{
								label: function(tooltipItem, data) {{
									return " "+values[tooltipItem.index]+" residues";
								}}
							  }}

						}}
						}}
					}});
				}}
			);



			// POINT MAP FOR 3PRED
d3.tsv("../prediction/3_prediction.csv").then(

				function(csv_data) {{
					var residues = csv_data.map(function(d) {{ return d.res; }});
					var flex = csv_data.map(function(d) {{ return d.pred_3; }});
					var pred = csv_data.map(function(d) {{ return d.P_max; }});
					var numberOfCharts = Math.ceil(residues.length / 50);
                    var idx = [];
                    for (var i=1; i <= residues.length; i++) {{
                        idx.push(i);
                    }}

					for (var i=1; i<=numberOfCharts; i++) {{
						var divC = document.createElement("div");
                        divC.classList.add("d-inline-flex");
						var canvas = document.createElement("canvas");
						canvas.id = "canvasChart3pred_"+i;
						canvas.height = 100;
						canvas.style = "margin-top: 25px";
						canvas.width = document.getElementById("principal_container").offsetWidth - 20;
						var mutate1_2 = false,
                            mutate3 = false,
                            mutate4 = false;
						if (residues.length > 50 && i == numberOfCharts && ((residues.length % 50)!=0)) {{
							if (residues.length % 50 >= 5) {{
								canvas.width = (document.getElementById("principal_container").offsetWidth - 30)*(residues.length % 50)*1.0 / 50 + 25*(1- (residues.length % 50)*1.0/50);
							}} else if (residues.length % 50 == 4) {{
                                    canvas.width = 135;
                                    canvas.style = "margin-left: -2.4rem !important";
                                    var mutate4 = true;
							}} else if (residues.length % 50 == 3) {{
								canvas.width = 135;
								canvas.style = "margin-left: -3.4rem !important";
								var mutate3 = true;
							}} else if (residues.length % 50 == 2 || residues.length % 50 == 1) {{
								canvas.width = 135;
								canvas.style = "margin-left: -4.6rem !important";
								var mutate1_2 = true;
							}}
						}}
						div = document.getElementById("divChart3pred");
						var arrow = document.createElement("img");
                        arrow.src = "../../../html/images/rigid_flexible_arrow.svg";
                        arrow.height = 60;
						arrow.style = "margin-left: -60px";
                        arrow.classList.add("align-self-center","mb-3");
						divC.appendChild(arrow);
                        divC.appendChild(canvas);
                        div.appendChild(divC);

						var chart = new Chart("canvasChart3pred_"+i, {{
							type: "line",
							data: {{
							  labels: residues.slice((i*50)-50, (i*50)),
							  datasets: [
											{{
												data: flex.slice((i*50)-50, (i*50)),
												flex: pred.slice((i*50)-50, (i*50)),
                                                idx: idx.slice((i*50)-50, (i*50)),
												tension: 0,
												fill: false,
												borderWidth: 0,
												pointBorderWidth: 0,
												pointRadius: 10,
												pointHoverRadius: 11,
												showLine: false,
												pointBackgroundColor: []
											}}
										]
							}},
							options: {{
								legend: {{
									display: false
								}},
								responsive: false,
								elements: {{
									point: {{
										pointStyle: "rect",
										borderWidth: 0,
										hoverBorderWidth: 0
									}}
								}},
								scales: {{
											xAxes: [{{
												ticks: {{
													autoSkip: false,
													maxRotation: 0,
													fontStyle: "bold",
													fontSize: 15,
													callback: function(value, index, values) {{
														if ((index+1) % 10 == 0 || (index+1) % 50 == 1) {{
															return [value,index+((i-1)*50)+1];
														}}
														else {{
															return value;
														}}
													}}
												}}
											}}],
											yAxes: [{{
												ticks: {{
													beginAtZero: true,
													fontSize: 16,
													stepSize: 1,
													autoSkip: true,
													padding: 10
												}}
											}}]
										}},
								tooltips: {{
									callbacks: {{
										title: function(tooltipItems, data) {{
                                            return data.labels[tooltipItems[0].index] + " (" + data.datasets[0].idx[tooltipItems[0].index] + ")";
                                        }},
                                        label: function(tooltipItem, data) {{
                                            return "Probabillity: "+((data.datasets[0].flex[tooltipItem.index]*100).toFixed(0))+" %";
                                        }}
									  }}

								}}
							}}
						}});

						if (mutate4){{
                            chart.options.layout.padding.left = 40;
                            chart.update();
                        }} else if (mutate3) {{
                            chart.options.layout.padding.left = 55;
                            chart.update();
                        }} else if (mutate1_2) {{
                            chart.options.layout.padding.left = 85;
                            chart.update();
                        }}

						var chartColors = {{
							color1: "rgba(35, 98, 118, 1)",
							color2: "rgba(0, 153, 159, 1)",
							color3: "rgba(0, 209, 166, 1)",
						}};


						var dataset = chart.data.datasets[0];
						for (var j = 0; j < dataset.data.length; j++) {{
							k = j + ((i-1)*50);
							if (pred[k] < 0.5) {{
								dataset.pointBackgroundColor[j] = chartColors.color3;
							}}
							else if ((pred[k] >= 0.5) && (pred[k] <= 0.6)){{
								dataset.pointBackgroundColor[j] = chartColors.color2;
							}}
							else {{
								dataset.pointBackgroundColor[j] = chartColors.color1;
							}}
						}}
						chart.update();

						window.chart = chart;
				}};
			}});





			// PIE CHART 5PRED
d3.tsv("../prediction/5_prediction.csv").then(

				function(csv_data) {{
					var flex = csv_data.map(function(d) {{ return d.pred_5; }});
					// parse class column, cast string into numerical value
					// and calculate occurences of each classes for pie chart
					var values = {{0: 0, 1: 0, 2: 0, 3: 0, 4: 0 }};
					for (var i=0; i < flex.length; i++) {{
						values[parseInt(flex[i], 10)] = (values[parseInt(flex[i], 10)] || 0) + 1;
					}};
					var theHelp = Chart.helpers;
					var chart_5pred = new Chart("pie_chart_5pred", {{
						type: "pie",
						data: {{
						  labels: ["0", "1", "2", "3", "4"],
						  datasets: [
										{{
											data: Object.values(values),
											backgroundColor: ["rgba(0, 78, 137, 1)", "rgba(35, 100, 170, 1)", "rgba(61, 165, 217, 1)", "rgba(254, 198, 1, 1)", "rgba(234, 115, 23, 1)"],
										}}
									]
						}},
						options: {{
                            legend: {{
                                display: true,
                                position: "right",
                                labels: {{
                                    generateLabels: function(chart) {{
                                        var data = chart.data;
                                        if (data.labels.length && data.datasets.length) {{
                                          return data.labels.map(function(label, i) {{
                                            var meta = chart.getDatasetMeta(0);
                                            var ds = data.datasets[0];
                                            var sum = ds["data"].reduce((pv, cv) => pv + cv, 0);
                                            var arc = meta.data[i];
                                            var custom = arc && arc.custom || {{}};
                                            var getValueAtIndexOrDefault = theHelp.getValueAtIndexOrDefault;
                                            var arcOpts = chart.options.elements.arc;
                                            var fill = custom.backgroundColor ? custom.backgroundColor : getValueAtIndexOrDefault(ds.backgroundColor, i, arcOpts.backgroundColor);
                                            var stroke = custom.borderColor ? custom.borderColor : getValueAtIndexOrDefault(ds.borderColor, i, arcOpts.borderColor);
                                            var bw = custom.borderWidth ? custom.borderWidth : getValueAtIndexOrDefault(ds.borderWidth, i, arcOpts.borderWidth);
                                              return {{
                                              // And finally :
                                              text: ((ds.data[i]*100)/sum).toFixed(0) + "% of class " + label,
                                              fillStyle: fill,
                                              strokeStyle: stroke,
                                              lineWidth: bw,
                                              hidden: isNaN(ds.data[i]) || meta.data[i].hidden,
                                              index: i
                                            }};
                                          }});
                                        }}
                                        return [];
                                    }}
                                }}
                            }},
                            title: {{
                                display: true,
                                fontFamily: "Roboto",
                                fontSize: 16,
                                fontColor: "black",
                                fontStyle: "normal",
                                lineHeight: 1.5,
								text: ["Global flexibility of the target sequence:", "residues are grouped into 5 predicted classes of flexibility"]
                            }},
							tooltips: {{
							callbacks: {{
								label: function(tooltipItem, data) {{
									return " "+values[tooltipItem.index]+" residues";
								}}
							  }}

						}}
						}}
					}});
				}}
			);

			// POINT MAP FOR 5PRED
d3.tsv("../prediction/5_prediction.csv").then(

				function(csv_data) {{
					var residues = csv_data.map(function(d) {{ return d.res; }});
					var flex = csv_data.map(function(d) {{ return d.pred_5; }});
					var pred = csv_data.map(function(d) {{ return d.P_max; }});
					var numberOfCharts = Math.ceil(residues.length / 50);
                    var idx = [];
                    for (var i=1; i <= residues.length; i++) {{
                        idx.push(i);
                    }}

					for (var i=1; i<=numberOfCharts; i++) {{
						var divC = document.createElement("div");
                        divC.classList.add("d-inline-flex");
						var canvas = document.createElement("canvas");
						canvas.id = "canvasChart5pred_"+i;
						canvas.height = 150;
						canvas.style = "margin-top: 25px";
						canvas.width = document.getElementById("principal_container").offsetWidth - 20;
						var mutate1_2 = false,
                            mutate3 = false,
                            mutate4 = false;
						if (residues.length > 50 && i == numberOfCharts && ((residues.length % 50)!=0)) {{
							if (residues.length % 50 >= 5) {{
								canvas.width = (document.getElementById("principal_container").offsetWidth - 30)*(residues.length % 50)*1.0 / 50 + 25*(1- (residues.length % 50)*1.0/50);
							}} else if (residues.length % 50 == 4) {{
                                    canvas.width = 135;
                                    canvas.style = "margin-left: -2.4rem !important";
                                    var mutate4 = true;
							}} else if (residues.length % 50 == 3) {{
								canvas.width = 135;
								canvas.style = "margin-left: -3.4rem !important";
								var mutate3 = true;
							}} else if (residues.length % 50 == 2 || residues.length % 50 == 1) {{
								canvas.width = 135;
								canvas.style = "margin-left: -4.6rem !important";
								var mutate1_2 = true;
							}}
						}}
						div = document.getElementById("divChart5pred");
						var arrow = document.createElement("img");
                        arrow.src = "../../../html/images/rigid_flexible_arrow.svg";
                        arrow.height = 100;
						arrow.style = "margin-left: -60px";
                        arrow.classList.add("align-self-center","mb-3");
						divC.appendChild(arrow);
                        divC.appendChild(canvas);
                        div.appendChild(divC);

						var chart = new Chart("canvasChart5pred_"+i, {{
							type: "line",
							data: {{
							  labels: residues.slice((i*50)-50, (i*50)),
							  datasets: [
											{{
												data: flex.slice((i*50)-50, (i*50)),
												flex: pred.slice((i*50)-50, (i*50)),
                                                idx: idx.slice((i*50)-50, (i*50)),
												tension: 0,
												fill: false,
												borderWidth: 0,
												pointBorderWidth: 0,
												pointRadius: 10,
												pointHoverRadius: 11,
												showLine: false,
												pointBackgroundColor: []
											}}
										]
							}},
							options: {{
								legend: {{
									display: false
								}},
								responsive: false,
								elements: {{
									point: {{
										pointStyle: "rect",
										borderWidth: 0,
										hoverBorderWidth: 0
									}}
								}},
								scales: {{
											xAxes: [{{
												ticks: {{
													autoSkip: false,
													maxRotation: 0,
													fontStyle: "bold",
													fontSize: 15,
													callback: function(value, index, values) {{
														if ((index+1) % 10 == 0 || (index+1) % 50 == 1) {{
															return [value,index+((i-1)*50)+1];
														}}
														else {{
															return value;
														}}
													}}
												}}
											}}],
											yAxes: [{{
												ticks: {{
													beginAtZero: true,
													fontSize: 16,
													stepSize: 1,
													autoSkip: true,
													padding: 10
												}}
											}}]
										}},
								tooltips: {{
									callbacks: {{
										title: function(tooltipItems, data) {{
                                            return data.labels[tooltipItems[0].index] + " (" + data.datasets[0].idx[tooltipItems[0].index] + ")";
                                        }},
                                        label: function(tooltipItem, data) {{
                                            return "Probabillity: "+((data.datasets[0].flex[tooltipItem.index]*100).toFixed(0))+" %";
                                        }}
									  }}

								}}
							}}
						}});

						if (mutate4){{
                            chart.options.layout.padding.left = 40;
                            chart.update();
                        }} else if (mutate3) {{
                            chart.options.layout.padding.left = 55;
                            chart.update();
                        }} else if (mutate1_2) {{
                            chart.options.layout.padding.left = 85;
                            chart.update();
                        }}


						var chartColors = {{
							color1: "rgba(35, 98, 118, 1)",
							color2: "rgba(0, 153, 159, 1)",
							color3: "rgba(0, 209, 166, 1)",
						}};


						var dataset = chart.data.datasets[0];
						for (var j = 0; j < dataset.data.length; j++) {{
							k = j + ((i-1)*50);
							if (pred[k] < 0.4) {{
								dataset.pointBackgroundColor[j] = chartColors.color3;
							}}
							else if ((pred[k] >= 0.4) && (pred[k] <= 0.5)){{
								dataset.pointBackgroundColor[j] = chartColors.color2;
							}}
							else {{
								dataset.pointBackgroundColor[j] = chartColors.color1;
							}}
						}}
						chart.update();

						window.chart = chart;
				}};
			}});







		</script>

		<!-- Footer -->
		<div id="footer"></div>
	</body></html>
    """)
