/* \filename init.js
   \author Yanxin Lu
   \date 04/04/2012
   
   This file is used for initialize all the html objects, including
   all the button events and etc.
*/

/* This function is used to initialize all the buttons in the page.
*/
function init_document() {
	$( document ).ready( function() {
		init_buttons();
	} );
}

/* Initialize all the buttons */
function init_buttons() {

	// cluster button
	$( "#run-button" ).click( function () {
		graph();
	} );

	// reset button
	$( "#reset-button" ).click( function () {
		reset();
	} );

	// download button
	$( "#download-button" ).click( function () {
		download();
	} );

	// button for showing downloading screen
	$( '#modal-button' ).click( function () {
		if( vis == null ) {
			show_message( "There is no network being shown" )
		} else {
			$( '#download-div' ).modal( 'show' );
		}

	} );
	$( '#modal-close-button' ).click( function () {
		$( '#download-div' ).modal( 'hide' );
	} );

	// we hide the output div first
	$( '#output-div' ).hide();

	// initialize the download div here
	$( '#download-div' ).modal();
	$( '#download-div' ).modal( 'hide' );

	// show a message first.
	show_message( "Please enter a query gene and select a species" );
}
