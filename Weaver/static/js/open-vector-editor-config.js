var editorProps = {
	shouldAutosave: false, //by default the editor does not autosave, setting this to true will trigger the onSave callback after any change to the sequenceData
	//supplying this function WILL make the editor FULLSCREEN BY DEFAULT
	handleFullscreenClose: () => {
		//do whatever you want here
		//UMD only:
		editor.close() //this calls reactDom.unmountComponent at the node you passed as the first arg
	},
	showMenuBar: true,
	showReadOnly: false, //default true
	disableSetReadOnly: true, //default false
	onSave: function(event, sequenceDataToSave, editorState, onSuccessCallback) {
		// console.info("event:", event);
		// console.info("sequenceData:", sequenceDataToSave);
		// console.info("editorState:", editorState);
		// To disable the save button after successful saving
		// either call the onSuccessCallback or return a successful promise :)
		// onSuccessCallback()
		//or
		return saveOVE(sequenceDataToSave);
	},
	onCopy: function(event, copiedSequenceData, editorState) {
		//the copiedSequenceData is the subset of the sequence that has been copied in the teselagen sequence format
		console.info("event:", event);
		console.info("sequenceData:", copiedSequenceData);
		console.info("editorState:", editorState);
		const clipboardData = event.clipboardData;
		clipboardData.setData("text/plain", copiedSequenceData.sequence);
		clipboardData.setData(
			"application/json",
			JSON.stringify(copiedSequenceData)
		);
		event.preventDefault();
		//in onPaste in your app you can do:
		// e.clipboardData.getData('application/json')
	},
	onPaste: function(event, editorState) {
		//the onPaste here must return sequenceData in the teselagen data format
		const clipboardData = event.clipboardData;
		let jsonData = clipboardData.getData("application/json")
		if (jsonData) {
			jsonData = JSON.parse(jsonData)
			if (jsonData.isJbeiSeq) {
				jsonData = convertJbeiToTeselagen(jsonData)
			}
		}
		const sequenceData = jsonData || {sequence: clipboardData.getData("text/plain")}
		return sequenceData
	},
	beforeAnnotationCreate: ({ //also works for edits (!)
		annotationTypePlural, //features/parts/primers
		annotation, //annotation info
		props //general props to the dialog
	}) => {
		//a handler to hook into when annotations (features/primers/parts) are being created
	},
	//regular click overrides, eg:
	featureClicked: ({annotation, event}) => {
		//do something here :)
	},
	// orf/primer/translation/cutsite/translationDouble/deletionLayer/replacementLayer/feature/part/searchLayer xxxxClicked can also be overridden

	rightClickOverrides: { //override what happens when a given feature/part/primer/translation/orf/cutsite/selectionLayer/lineageLine gets right clicked
		//the general format is xxxxRightClicked eg:
		selectionLayerRightClicked: (items, {annotation}, props) => {
			return [...items, {
				//props here get passed directly to blueprintjs MenuItems
				text: "Create Part",
				onClick: () => console.info('hey!≈')
			}]
		}
	},
	PropertiesProps: {
		// the list of tabs shown in the Properties panel
		propertiesList: [
			"general",
			"features",
			"parts",
			"primers",
			"translations",
			"cutsites",
			"orfs",
			"genbank"
		]
	},
	ToolBarProps: {
			toolList: [
				"saveTool",
				//you can override a tool like so:
				{name: "downloadTool", Dropdown: () => {
					return "Pronto!"
				}},
				"importTool",
				"undoTool",
				"redoTool",
				"cutsiteTool",
				"featureTool",
				"alignmentTool",
				// "oligoTool",
				"orfTool",
				// "viewTool",
				"editTool",
				"findTool",
				"visibilityTool"
				// "propertiesTool"
			]
		},
		StatusBarProps: {
			//these are the defaults:
			showCircularity: true,
			showReadOnly: false,
			showAvailability: false
		},
		onDigestSave: () => {} //tnr: NOT YET IMPLEMENTED
};
var editorState = {
	//note, sequence data passed here will be coerced to fit the Teselagen data model (Teselagen JSON)
	sequenceData: {
	//Open Vector Editor data model
		sequence: "atagatagagaggcccg",
		features: [
			{
				color: "#b3b3b3", //you can override the default color for each individual feature if you want
				type: "misc_feature",
				start: 0, //start and end are 0-based inclusive for all annotations
				end: 10,
				id: 'yourUniqueID',
				forward: true //ie true=positive strand     false=negative strange
			}
		],
		parts: []
	},
	sequenceDataHistory: {}, //clear the sequenceDataHistory if there is any left over from a previous sequence
	annotationVisibility: {
		features: true
	},
	panelsShown: [
		[
			{
				id: "sequence",
				name: "Sequence Map",
				active: true
			}
		],
		[
			{
				id: "circular",
				name: "Circular Map",
				active: true
			},
			{
				id: "rail",
				name: "Linear Map",
				active: false
			},
			{
				id: "properties",
				name: "Properties",
				active: false
			}
		]
	],
};

async function sequenceToJson() {
    if(typeof sequence_file_contents == 'undefined') {return;}
    var jsonOutput = await window.bioParsers.anyToJson(sequence_file_contents);
    if(jsonOutput[0]['success']){
        editorState['sequenceData'] = jsonOutput[0]['parsedSequence'];
        editor = window.createVectorEditor(document.getElementById("ove-viewer"), editorProps);
        editor.updateEditor(editorState);
    } else {
        window.toastr.success("Error parsing plasmid sequence file.");
    }
}
sequenceToJson();