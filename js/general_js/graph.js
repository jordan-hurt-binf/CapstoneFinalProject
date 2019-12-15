/* \filename complexfinder_cluster.js
   \ written  by Yanxin Lu  		11/12/2011  
   \ modified by Ryan Atanasoff 	11/18/2017
   \ modified by London Steele, James Yen, Lily Wise, Abby Antrich Fall 2018
   \ modified by Group 2 Fall 2019
   This script provides the main functions for the SSFinder home page.
*/

// global variable for visualization object
var vis;

// used to expand divs
var add = 0;

// a variable to check whether we are showing tabs
var showing_tabs = false;

// default height for bottom div
var default_height = 530;

/* This function gathers all the information from the pages and send
 * them to the server and let the server to do the graphing. It then
 * receives the result and displays it. This should be the main
 * function in this file since the web calls this function directly.
 */
function graph() {
	// get all the graph options

	// get both of the genes
	var gene_1 = $('#gene1').val().toUpperCase();
	var gene_2 = $('#gene2').val().toUpperCase();

	// checking that values were given to the variables
	 if( gene_1 == "" || gene_2 == "" ) {
	 	// error message if both values were not initialized
	 	show_message( "Invalid Input: Please enter two genes." );
	 	return;
	 } 


	// get the network from the server
	send_cluster_request(gene_1, gene_2);
}

/* This function shows a message on the page */
function show_message( msg ){
        show_message(msg, false);
}

function show_message( msg, error) {
	if( msg == "Success!" ) {
		$( '#message-div' ).attr( 'class', 'alert alert-success' );
	} else if (error == true){
		$( '#message-div' ).attr( 'class', 'alert alert-error' );

	} else {
		$( '#message-div' ).attr( 'class', 'alert' );
	}
	// show the message in the description div
	$( "#message-div" ).html( msg );
}

/* This function takes a network in json format and draws it on the web */
function draw( network, gene_1, gene_2 ) {

	// hide the description div first
	$( '#description-div' ).hide();
	// show the success message
	show_message( "Success!" );

	// set up the output-div for displaying
	if( !$( '#output-div' ).is( ":visible" ) ) {
		$( '#output-div' ).fadeIn( "fast" );
	}

	$("#info-div").show();

	// clear the tables and tabs
	$( '#tab1' ).html( '' );
	$( '#tab2' ).html( '' );
	$( '#tab-content1' ).html( '' );
	$( '#tab-content2' ).html( '' );
	if( showing_tabs ) {
		showing_tabs = false;
		remove_tabs();
	}

	
	// init and draw
	vis = cytoscape({
		container: document.getElementById('network-div'),
		style: cytoscape.stylesheet()
			.selector('node[type = "foo"]')
			.css({
				'label': 'data(id)',
				'text-valign': 'center',
				'shape': 'circle',
				'width': '70px',
				'height': '70px',
				'background-color': '#DFF0D8',
				'color': '#468847',
				'border-width': '1px',
				'border-color': '#468847',
				'font-size': 14,
				'font': '300 14px "Helvetica Neue", Helvetica'
			})
			.selector('node[type = "noFoo"]')
			.css({
				'label': 'data(id)',
				'text-valign': 'center',
				'shape': 'circle',
				'width': '70px',
				'height': '70px',
				'background-color': ' #F7FCFF',
				'color': '#0088CC',
				'border-width': '1px',
				'border-color': '#0088CC',
				'font-size': 14,
				'font': '300 14px "Helvetica Neue", Helvetica'
			})
			.selector('node.highlighted')
			.css({
				'border-color': '#2E64FE'
			})
			.selector('edge.highlighted')
			.css({
				'line-color': '#2E64FE'
			})
			.selector('edge')
				.css({
				'width': 1.5,
				'curve-style': 'haystack'
			}),
	});
	/*display options for vis object*/
	var options = {
		name: 'concentric',
		
		fit: true, // whether to fit the viewport to the graph
		padding: 30, // the padding on fit
		startAngle: 3 / 2 * Math.PI, // where nodes start in radians
		clockwise: true, // whether the layout should go clockwise (true) or counterclockwise/anticlockwise (false)
		equidistant: false, // whether levels have an equal radial distance betwen them, may cause bounding box overflow
		minNodeSpacing: 20, // min spacing between outside of nodes (used for radius adjustment)
		boundingBox: undefined, // constrain layout bounds; { x1, y1, x2, y2 } or { x1, y1, w, h }
		avoidOverlap: true, // prevents node overlap, may overflow boundingBox if not enough space
		nodeDimensionsIncludeLabels: false, // Excludes the label when calculating node bounding boxes for the layout algorithm
		height: undefined, // height of layout area (overrides container height)
		width: undefined, // width of layout area (overrides container width)
		spacingFactor: undefined, // Applies a multiplicative factor (>0) to expand or compress the overall area that the nodes take up
		concentric: function( node ){ // returns numeric value for each node, placing higher nodes in levels towards the centre
		return node.degree();
		},
		levelWidth: function( nodes ){ // the letiation of concentric values in each level
		return nodes.maxDegree() / 4;
		},
		ready: undefined, // callback on layoutready
		stop: undefined, // callback on layoutstop 
	};

	/* sets check equal to official or systematic name of what the user inputed  */
	//var check = gene_1;
	//var cur_off = "";
	//for (var key in var key in network[0]) {
	    /* CHECKS IF THE VALUE IS A SYNONYM OR NOT AND MAY RESET THE VALUE */
	//}

	/*add nodes to vis object*/	
	for (var key in network.nodes) {

		if( key == gene_1 || key == gene_2 ) {
			vis.add({
				data: { 
					id: key,
					type: "foo",
					process: network.nodes[key].process,
					synonyms: network.nodes[key].synonyms
				}
			});
		}
		else {
			vis.add({
				data: { 
					id: key,
                     			type: "noFoo",
                    			process: network.nodes[key].process,
                   			synonyms: network.nodes[key].synonyms
				}
			});
		}
	}

	/*add edges to vis object*/
	for (var key in network.edges){
		vis.add({
			data: {
				id: key,
				source: network.edges[key][0],
				target: network.edges[key][1],
				ss: network.edges[key][2]
			}
		});
	}
	
	/*draw graph, check zoom to see if divs need expanding*/
	vis.layout( options );
	var zoom_temp = vis.zoom();
	vis.panzoom({});
	vis.zoom(.9);
	while(zoom_temp < .9){
		zoom_temp+=.1;
		add += 50;
	}
	vis.userZoomingEnabled( true );

	/*set initial display to user's gene*/
	for(var key in vis._private.elements._private.ids){
		
		var temp = vis._private.elements._private.ids[key]._private.data;
		print_gene_info(temp);
		if( !showing_tabs ) {
			construct_tabs();
			showing_tabs = true;
		}
		break;
		//}
	}

	/*on click effects*/
	var previousEvent;
	vis.on('click', 'node', function(e){
		if( e.cyTarget._private.group == "nodes" ) {
			if(previousEvent !== undefined){
				previousEvent.removeClass('highlighted');
			}
			print_gene_info( e.cyTarget._private.data );
			e.cyTarget.addClass('highlighted');
			previousEvent = e.cyTarget;
			if( !showing_tabs ) {
				construct_tabs();
				showing_tabs = true;
			}
	}});
	vis.on('click', 'edge', function(e){
		if( e.cyTarget._private.group = "edges") {
			if(previousEvent !== undefined){
				previousEvent.removeClass('highlighted');
			}
			console.log(e.cyTarget._private.data);
			print_edge_info( e.cyTarget._private.data );
			e.cyTarget.addClass('highlighted');
			previousEvent = e.cyTarget;
			if( !showing_tabs ) {
				construct_tabs();
				showing_tabs = true;
			}
		}
	});
	temp_function2(add);
}

/* This function creates a ajax request */
function send_cluster_request( gene_1, gene_2 ) {
	var count=60;

	var counter;
	function timer() {
		count=count-1;
		if (count <= 0) {
			clearInterval(counter);
                        return;
		}
		document.getElementById("modalcd").innerHTML='Response no later than in '+ count+' second.';
		show_message('Discovering the functional pathways between the query genes. Response no later than in '+ count + ' seconds.');
	}
	document.getElementById("modalcd").innerHTML='Response no later than in '+ count+' second.';
	show_message('Discovering the functional pathways between the query genes. Response no later than in '+ count + ' seconds.');

	counter=setInterval(timer, 1000); //1000 will  run it every 1 second
	$('#countdownmodal').modal();
	// send a post request to the server
	let url = '/python/test.py';
	$.ajax({ url, // url of the script
		// all the data we want to send
		data: { 
			startGene: gene_1,
			endGene: gene_2,
		},
		timeout: 60000, // timeout
		dataType: "json", // datatype
		type: "POST", // this is a post request

		// call back functions
		success: function( network ) {
			console.log(network);
			clearInterval(counter);
			$('#countdownmodal').modal('hide');
			if( network.nodes.length > 100 ) {
				show_message( "Too many nodes. Cannot be visualized." , true);
				$('#countdownmodal').modal('hide');
				clear();
			} else {
				// draw the graph
				$('#countdownmodal').modal('hide');
				draw( network, gene_1, gene_2);
			}},
		error: function( jqXHR, error ) {
			if( error == 'Your query timed out' ) {
				clearInterval(counter);
				$('#countdownmodal').modal('hide');
				show_message( "Timeout.", true );
			}else{
				console.log(jqXHR,error);
				clearInterval(counter);
				$('#countdownmodal').modal('hide');
				show_message("An unexpected error (" + error + ") occured. jqXHR: " + jqXHR.responseText + " Please try again.", true);
			}
		},
		async: false
	});
}

// hide the output and info division
function clear(){
	$( "#output-div").hide();
	$( "#info-div").hide(); 
}

/* This function resets everything. */
function reset() {

	// hide all other divs
	$( '#output-div' ).hide();
	show_message( "Please enter two query genes." );

	$( "#gene" ).val("");
	$( "#param" ).val("20");

	$( '#bottom-div' ).height( default_height + 70);

	// show the description-div
	$( '#description-div' ).fadeIn( "fast" );

	// nullify the visualization object
	vis = null;
}

/* This function returns the text for description */
function get_description() {
	return "Description Goes here";
}

/* This function handles the downloading function of the website. It
   first checks whether the visualization object is there, and if
   it is, it will call the exportNetwork function to send the image
   data to the php script we have on the server in order to le the
   users to download the network image */
function download() {
	if( vis != null ) {
		var format = $( '#format' ).val();
		// this is for downloading
		if( format == "png" ) {
			var png64 = vis.png();
			$('#png-eg').attr('src', png64);
		} else if( format == "svg" ) {
			vis.exportNetwork( "svg", "/cgi-bin/complexfinder/export.php?type=svg" );
		} else if( format == "pdf" ) {
			vis.exportNetwork( "pdf", "/cgi-bin/complexfinder/export.php?type=pdf" );
		} else if( format == "xgmml" || format == "graphml" ) {
			vis.exportNetwork( format, "/cgi-bin/complexfinder/export.php?type=xml" );
		} else if( format == "sif" ) {
			vis.exportNetwork( format, "/cgi-bin/complexfinder/export.php?type=txt" );
		}
	}
}

/* This function modify the network_info's html content to display the gene
   data by constructing tables. */
function print_gene_info( data ) {
	$( '#tab1' ).html( '' );
	$( '#tab2' ).html( '' );
	$( '#tab-content1' ).html( '' );
	$( '#tab-content2' ).html( '' );

	content = "";

	content += "<table class=\"table\">\n";

	content += "<tr><th colspan=\"3\">Synonyms of " + data.id + "</th></tr>\n";
	content += "</tr>\n";

	content += "<tbody>\n" 
	content += "<tr>\n";

	if (data.synonyms.length > 0){
		// Printing the synonyms for the gene
		for (i = 0; i < data.synonyms.length; i++){
			content += "<td>";
			content += data.synonyms[i];
			content += "</td>";
			content += "<tr>\n";
		}

	} else {
		content += "<td>" + data.id + " has no synonyms</td>";
	}
		

	content += "</tr>\n";
	content += "</tbody>\n";
	
	// Printing the Links to further Information
	content += "<tr><th colspan=\"3\">Links to More Information</th></tr>\n";
	content += "</tr>\n";

	content += "<tbody>\n" 

	content += "<tr>\n<td>";
	content += "<a href='https://www.ncbi.nlm.nih.gov/gene/?term=" + data.id + "' target='_blank'>Entrez Gene Query Results</a>";
	content += "</td>";
	content += "</tr>\n"; 

	content += "<tr>\n<td>";
	content += "<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=" + data.id + "&keywords=" + data.id + "' target='_blank'>GeneCards</a>";
	content += "</td>";
	content += "</tr>\n";

	content += "</tbody>\n";
	content += "</table>\n";

	// SOME FORMATTING EXAMPLE
    	$( '#tab1' ).html( "General Gene Information" );
	$( '#tab-content1' ).html( content );
   
    
    	// adjust the height of bottom-div
    	$( '#network-div' ).height( default_height);
    	$( '#bottom-div' ).height( Math.max( default_height, $( '#info-tabs' ).height() ) + 70 );
    	$( '#info-div' ).height( Math.max( default_height, $( '#info-tabs' ).height() ) );
}

function temp_function(){
	$( '#network-div' ).height( default_height );
	$( '#bottom-div' ).height( Math.max( default_height, $( '#tab-content2' ).height() ) + 134 );
	$( '#info-div' ).height( Math.max( default_height, $( '#tab-content2' ).height() ) + 64);
}

function temp_function1(){
        $( '#network-div' ).height( default_height );
        var h = $( '#tab-content1' ).height() + 92;
        $( '#info-div' ).height( Math.max( h, default_height ) );
        $( '#bottom-div' ).height( $( '#info-div' ).height() + 70 );

}

function temp_function2(add) {
	if($( '#network-div' ).height()+add < 550 ){
		$( '#network-div' ).height( default_height + add );
	}
	else{
		$( '#network-div' ).height( 550 );
	}
	if($( '#network-div' ).width()+add+40 < 600 ){
		$( '#network-div' ).width( default_height + add );
	}
	else{
		$( '#network-div' ).width( 580 );
		if(0 < add){
			vis.fit();
		}
		else{
			vis.zoom(.9);
			vis.center();
		}
	}
	$( '#bottom-div' ).height( $( '#network-div' ).height() + 70);
	default_height = $( '#bottom-div' ).height() - 70;
}

/* This function prints the edge info */
function print_edge_info( data ) {
	$( '#tab1' ).html( '' );
	$( '#tab2' ).html( '' );
	$( '#tab-content1' ).html( '' );
	$( '#tab-content2' ).html( '' );

	//SOME FORMATTING EXAMPLE
	$( '#tab1' ).html( "General Edge Information" );

	var content = "";
	// print the general info
	// ----------------------

	// create the table header
	content += "<table class=\"table\">\n";

	content += "<tr><th colspan=\"3\">Edge between " + data.source + " and " + data.target + "</th></tr>\n";
	content += "</tr>\n";

	// create the body
	content += "<tbody>\n";
	content += "<tr><td>Source</td><td>" + data.source  + "</td></tr>\n";
	content += "<tr><td>Target</td><td>" + data.target + "</td<</tr>\n";
	content += "</tbody>\n";

	content += "<tbody>\n";
	content += "<tr><td>Similarity Score</td><td>" + data.ss  + "</td></tr>\n";
	content += "</tbody>\n";
	
	// finish this table
	content += "</table>\n";

	$( '#tab-content1' ).html( content );


	// adjust the height of bottom-div
	$( '#network-div' ).height( default_height);
	$( '#bottom-div' ).height( Math.max( default_height, $( '#info-tabs' ).height() ) + 70 );
	$( '#info-div' ).height( Math.max( default_height, $( '#info-tabs' ).height() ) );
}

function construct_tabs() {
	// add info-tabs classes so that it actually uses the Bootstrap UI
	$( '#info-tabs' ).addClass( "tabbable tabs-below" );

	// make the actual tabs
	$( '#info-tabs ul:first-child' ).addClass( "nav nav-tabs" );
}

function remove_tabs() {
	$( '#info-tabs' ).removeClass( "tabbable tabs-below" );
	$( '#info-tabs ul:first-child' ).removeClass( "nav nav-tabs" );
}
