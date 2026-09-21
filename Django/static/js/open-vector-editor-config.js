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
	disableSetReadOnly: false, //default false
	onSave: function(event, sequenceDataToSave, editorState, onSuccessCallback) {
		// console.info("event:", event);
		// console.info("sequenceData:", sequenceDataToSave);
		// console.info("editorState:", editorState);
		// To disable the save button after successful saving
		// either call the onSuccessCallback or return a successful promise :)
		// onSuccessCallback()
		//or
		return saveOVE(weaverSequenceDataForSave(sequenceDataToSave));
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
	onSelectionOrCaretChanged: function(selectionState) {
		const selection = selectionState ? selectionState.selectionLayer : null;
		if (selection && selection.start > -1 && selection.end > -1) {
			weaverLatestSelection = {
				start: Number(selection.start),
				end: Number(selection.end),
				overlapsSelf: Boolean(selection.overlapsSelf || selection.start > selection.end)
			};
			editorState.selectionLayer = {...selection};
		} else {
			weaverLatestSelection = null;
		}
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
			const selection = annotation || (props ? props.selectionLayer : null);
			const hasSelection = selection && selection.start > -1 && selection.end > -1;
			return [...items, {
				text: "Design PCR",
				disabled: !hasSelection || typeof plasmid_pcr_design_path === "undefined",
				onClick: () => {
					if (!hasSelection || typeof plasmid_pcr_design_path === "undefined") {
						return;
					}
					const params = new URLSearchParams({
						start: selection.start,
						end: selection.end,
						margin: 300
					});
					window.location.href = plasmid_pcr_design_path + "?" + params.toString();
				}
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
var weaverAmpliconCandidates = [];
var weaverNonOverlappingAmplicons = [];
var weaverAmpliconCandidatesLoaded = false;
var weaverVisibleAmpliconIds = new Set();
var weaverExpandedDimerIds = new Set();
var weaverSelectedPolymeraseKey = "sapphire_amp";
var weaverDigestResults = [];
var weaverVisibleDigestId = "";
var weaverLatestSelection = null;
var weaverDigestOriginalSequenceData = null;
var weaverDigestRotation = 0;
var weaverAmpliconOriginalSequenceData = null;
var weaverAmpliconRotation = 0;

var WEAVER_POLYMERASES = {
	sapphire_amp: {
		name: "SapphireAmp",
		minSecondsUnderBp: 1000,
		minSeconds: 5,
		secondsPerKb: 10,
		roundUpToSeconds: 5
	},
	taq: {
		name: "Taq",
		minSecondsUnderBp: 0,
		minSeconds: 0,
		secondsPerKb: 60,
		roundUpToSeconds: 5
	},
	phusion: {
		name: "Phusion",
		minSecondsUnderBp: 0,
		minSeconds: 0,
		secondsPerKb: 30,
		roundUpToSeconds: 5
	},
	custom_polymerase: {
		name: "Custom polymerase",
		minSecondsUnderBp: 0,
		minSeconds: 0,
		secondsPerKb: null,
		roundUpToSeconds: 5
	}
};
const WEAVER_NO_POLYMERASE_OPTION = "__no_polymerase__";
const WEAVER_CUSTOM_POLYMERASE_KEY = "custom_polymerase";
const WEAVER_ANNEALING_TEMPERATURE_HELP = "Recommended annealing temperature: Ta Opt = 0.3 x Tm of the less stable primer + 0.7 x Tm of the PCR product - 14.9.";

function cloneForWeaver(value) {
	return JSON.parse(JSON.stringify(value || {}));
}

function setPrimerMatchStatus(message, isError) {
	const status = document.getElementById("ove-primer-match-status");
	if (!status) {
		const panelSummary = document.getElementById("ove-amplicon-panel-summary");
		if (panelSummary && message && isError) {
			panelSummary.textContent = message;
		}
		return;
	}
	status.textContent = message || "";
	status.classList.toggle("text-danger", Boolean(isError));
	status.classList.toggle("text-secondary", !isError);
}

function setWeaverSidePanelOpen(panelName) {
	const ampliconPanel = document.getElementById("ove-amplicon-panel");
	const digestPanel = document.getElementById("ove-digest-panel");
	if (!ampliconPanel || !digestPanel) {
		return;
	}
	ampliconPanel.classList.toggle("is-collapsed", panelName !== "amplicons");
	digestPanel.classList.toggle("is-collapsed", panelName !== "digests");
	ampliconPanel.classList.toggle("is-peer-open", panelName === "digests");
	digestPanel.classList.toggle("is-peer-open", panelName === "amplicons");
}

function setAmpliconPanelOpen(isOpen) {
	setWeaverSidePanelOpen(isOpen ? "amplicons" : "");
}

function setDigestPanelOpen(isOpen) {
	setWeaverSidePanelOpen(isOpen ? "digests" : "");
}

function stopCircularMapWheelPageScroll() {
	const viewer = document.getElementById("ove-viewer");
	if (!viewer || viewer.dataset.weaverCircularWheelGuard === "true") {
		return;
	}
	viewer.dataset.weaverCircularWheelGuard = "true";
	viewer.addEventListener("wheel", function(event) {
		if (!event.target || !event.target.closest(".veCircularView")) {
			return;
		}
		event.preventDefault();
		event.stopPropagation();
	}, {passive: false});
}

	function valueFromNotes(annotation, key) {
		const value = annotation && annotation.notes ? annotation.notes[key] : null;
		if (Array.isArray(value)) {
			return value[0] || "";
		}
		return value || "";
	}

	function valuesFromNotes(annotation, key) {
		const value = annotation && annotation.notes ? annotation.notes[key] : null;
		if (Array.isArray(value)) {
			return value.filter((item) => String(item || "").trim());
		}
		return value ? [value] : [];
	}

function ampliconCoordinateLabel(annotation) {
	if (!annotation) {
		return "";
	}
	const start = Number(annotation.start) + 1;
	const end = Number(annotation.end) + 1;
	return annotation.overlapsSelf ? start + ".." + end + " circular" : start + ".." + end;
}

function primerIdLabel(annotation) {
	const fwdId = valueFromNotes(annotation, "fwd_primer_id");
	const revId = valueFromNotes(annotation, "rev_primer_id");
	if (!fwdId && !revId) {
		return "";
	}
	return "IDs: " + (fwdId || "?") + " + " + (revId || "?");
}

function annealingTemperatureLabel(annotation) {
	const annealingTm = valueFromNotes(annotation, "recommended_annealing_tm");
	return annealingTm ? "Ta " + annealingTm + " °C" : "";
}

	function annealingTemperatureTooltip(annotation) {
		const fwdTm = valueFromNotes(annotation, "fwd_tm");
		const revTm = valueFromNotes(annotation, "rev_tm");
		const productTm = valueFromNotes(annotation, "product_tm");
	const details = [];
	if (fwdTm) {
		details.push("FWD Tm: " + fwdTm + " °C");
	}
	if (revTm) {
		details.push("REV Tm: " + revTm + " °C");
	}
		if (productTm) {
			details.push("Product Tm: " + productTm + " °C");
		}
		return details.length ? WEAVER_ANNEALING_TEMPERATURE_HELP + "\n" + details.join("\n") : WEAVER_ANNEALING_TEMPERATURE_HELP;
	}

		function primerComplementarityLabel(annotation) {
			const primer3Risk = valueFromNotes(annotation, "primer3_dimer_risk");
			if (primer3Risk) {
				return "dimer: " + primer3Risk.toLowerCase();
			}
			const severity = valueFromNotes(annotation, "primer_complementarity_severity");
			if (!severity) {
				return "";
			}
			return "dimer: " + severity;
		}

		function primerComplementarityWarnings(annotation) {
			return [
				...valuesFromNotes(annotation, "primer_complementarity_warnings"),
				valueFromNotes(annotation, "primer3_dimer_reasons"),
				valueFromNotes(annotation, "primer3_dimer_recommendation"),
				valueFromNotes(annotation, "primer3_dimer_error")
			].filter((warning) => String(warning || "").trim());
		}

	function primerComplementarityAlignment(annotation) {
		const lines = [
			valueFromNotes(annotation, "primer_complementarity_alignment_fwd"),
			valueFromNotes(annotation, "primer_complementarity_alignment_match"),
			valueFromNotes(annotation, "primer_complementarity_alignment_rev")
		].filter((line) => String(line || "").trim());
		return lines.length === 3 ? lines : [];
	}

		function primerComplementaritySummary(annotation) {
			const primer3Risk = valueFromNotes(annotation, "primer3_dimer_risk");
			if (primer3Risk) {
				const dimerTm = valueFromNotes(annotation, "primer3_dimer_tm");
				const dimerDg = valueFromNotes(annotation, "primer3_dimer_dg");
				const longestRun = valueFromNotes(annotation, "primer3_dimer_longest_run") || "0";
				const threePrimeRun = valueFromNotes(annotation, "primer3_dimer_3prime_run") || "0";
				const details = [
					"PCR dimer risk " + primer3Risk + ".",
					dimerTm ? "Dimer Tm " + dimerTm + " °C." : "",
					dimerDg ? "ΔG " + dimerDg + " kcal/mol." : "",
					"Longest run " + longestRun + " bp; 3' run " + threePrimeRun + " bp."
				].filter((detail) => detail);
				return details.join(" ");
			}
			const severity = valueFromNotes(annotation, "primer_complementarity_severity") || "none";
			const maxMatch = valueFromNotes(annotation, "primer_complementarity_max") || "0";
			const threePrime = valueFromNotes(annotation, "primer_complementarity_3prime") || "0";
			const bothThreePrime = valueFromNotes(annotation, "primer_complementarity_both_3prime") || "0";
			return "Severity " + severity + ". Longest match " + maxMatch + " bp; 3' match " + threePrime + " bp; both 3' ends " + bothThreePrime + " bp.";
		}

		function primer3DimerStructure(annotation) {
			return valueFromNotes(annotation, "primer3_dimer_structure");
		}

	function togglePrimerComplementarityDetails(annotation) {
		if (!annotation || !annotation.id) {
			return;
		}
		if (weaverExpandedDimerIds.has(annotation.id)) {
			weaverExpandedDimerIds.delete(annotation.id);
		} else {
			weaverExpandedDimerIds.add(annotation.id);
		}
		renderAmpliconPanel(weaverAmpliconCandidates);
	}

function rangeForAmpliconSelection(annotation) {
	const start = Number(annotation.start);
	const end = Number(annotation.end);
	const sequenceLength = editorState.sequenceData && editorState.sequenceData.sequence ?
		editorState.sequenceData.sequence.length :
		0;
	if ((annotation.overlapsSelf || end < start) && sequenceLength) {
		return [
			{
				start,
				end: sequenceLength - 1,
				overlapsSelf: false
			},
			{
				start: 0,
				end,
				overlapsSelf: false
			}
		];
	}
	return {
		start,
		end,
		overlapsSelf: false
	};
}

function selectedPolymerase() {
	return WEAVER_POLYMERASES[weaverSelectedPolymeraseKey] || null;
}

function roundUpToStep(value, step) {
	const safeStep = Number(step) > 0 ? Number(step) : 1;
	return Math.ceil(Number(value) / safeStep) * safeStep;
}

function productSizeFromAmplicon(amplicon) {
	return Number(valueFromNotes(amplicon, "product_size")) || 0;
}

function extensionSecondsForAmplicon(amplicon) {
	const polymerase = selectedPolymerase();
	const productSize = productSizeFromAmplicon(amplicon);
	if (!polymerase || !productSize) {
		return null;
	}
	if (productSize < polymerase.minSecondsUnderBp) {
		return polymerase.minSeconds;
	}
	const rawSeconds = productSize / 1000 * polymerase.secondsPerKb;
	return roundUpToStep(rawSeconds, polymerase.roundUpToSeconds);
}

	function extensionTimeLabel(amplicon) {
	const seconds = extensionSecondsForAmplicon(amplicon);
	if (seconds === null) {
		return "";
	}
	if (seconds < 60) {
		return "Ext " + seconds + " s";
	}
	const minutes = Math.floor(seconds / 60);
	const remainder = seconds % 60;
		return remainder ? "Ext " + minutes + " min " + remainder + " s" : "Ext " + minutes + " min";
	}

	function partialPrimerBindingLabel(amplicon) {
		return valueFromNotes(amplicon, "partial_primer_binding") === "true" ? "3' partial" : "";
	}

	function partialPrimerBindingTooltip(amplicon) {
		const details = [];
		if (valueFromNotes(amplicon, "fwd_partial_binding") === "true") {
			details.push(
				"FWD: " +
				valueFromNotes(amplicon, "fwd_binding_length") +
				" bp aligned; 5' unaligned " +
				valueFromNotes(amplicon, "fwd_unmatched_5")
			);
		}
		if (valueFromNotes(amplicon, "rev_partial_binding") === "true") {
			details.push(
				"REV: " +
				valueFromNotes(amplicon, "rev_binding_length") +
				" bp aligned; 5' unaligned " +
				valueFromNotes(amplicon, "rev_unmatched_5")
			);
		}
		return details.length ? details.join("\n") : "Primer binding uses the 3' end only.";
	}

	function fallbackCopyText(text) {
		const textarea = document.createElement("textarea");
		textarea.value = text;
		textarea.setAttribute("readonly", "");
		textarea.style.position = "fixed";
		textarea.style.left = "-9999px";
		document.body.appendChild(textarea);
		textarea.select();
		try {
			document.execCommand("copy");
			return true;
		} finally {
			document.body.removeChild(textarea);
		}
	}

	async function copyAmpliconSequence(amplicon) {
		const sequence = valueFromNotes(amplicon, "amplicon_sequence");
		if (!sequence) {
			setPrimerMatchStatus("Amplicon sequence is not available.", true);
			return;
		}
		try {
			if (navigator.clipboard && navigator.clipboard.writeText) {
				await navigator.clipboard.writeText(sequence);
			} else if (!fallbackCopyText(sequence)) {
				throw new Error("Clipboard fallback failed");
			}
			setPrimerMatchStatus("Amplicon copied to clipboard.");
		} catch (error) {
			setPrimerMatchStatus("Unable to copy amplicon.", true);
		}
	}

	function blastAmpliconSequence(amplicon) {
		if (typeof services_blast_path === "undefined") {
			setPrimerMatchStatus("BLAST service is not available.", true);
			return;
		}
		const sequence = valueFromNotes(amplicon, "amplicon_sequence");
		if (!sequence) {
			setPrimerMatchStatus("Amplicon sequence is not available.", true);
			return;
		}
		const forwardPrimer = valueFromNotes(amplicon, "fwd_primer") || valueFromNotes(amplicon, "fwd_primer_id");
		const reversePrimer = valueFromNotes(amplicon, "rev_primer") || valueFromNotes(amplicon, "rev_primer_id");
		const blastName = forwardPrimer && reversePrimer ?
			forwardPrimer + " + " + reversePrimer :
			(amplicon.name || forwardPrimer || reversePrimer || "Amplicon");
		const params = new URLSearchParams({
			name: blastName,
			sequence,
			project: "a"
		});
		window.open(services_blast_path + "?" + params.toString(), "_blank", "noopener");
	}

function loadCustomPolymeraseTime() {
	try {
		const secondsPerKb = Number(localStorage.getItem("weaverCustomPolymeraseSecondsPerKb"));
		if (secondsPerKb > 0) {
			WEAVER_POLYMERASES[WEAVER_CUSTOM_POLYMERASE_KEY].secondsPerKb = secondsPerKb;
		}
	} catch (error) {
	}
}

function firstPolymeraseKey() {
	return Object.keys(WEAVER_POLYMERASES)[0] || "";
}

function promptForPolymeraseSecondsPerKb(key) {
	const polymerase = WEAVER_POLYMERASES[key] || null;
	if (!key || !polymerase) {
		return false;
	}
	const secondsPerKbInput = window.prompt(
		"Extension seconds per kb",
		polymerase.secondsPerKb ? String(polymerase.secondsPerKb) : ""
	);
	const secondsPerKb = Number(secondsPerKbInput);
	if (!secondsPerKbInput) {
		return false;
	}
	if (!secondsPerKb || secondsPerKb <= 0) {
		window.alert("Extension seconds per kb must be a positive number.");
		return false;
	}
	WEAVER_POLYMERASES[key] = {...polymerase, secondsPerKb};
	if (key === WEAVER_CUSTOM_POLYMERASE_KEY) {
		localStorage.setItem("weaverCustomPolymeraseSecondsPerKb", String(secondsPerKb));
	}
	return true;
}

function populatePolymeraseSelect(selectElement) {
	selectElement.innerHTML = "";
	const polymeraseEntries = Object.entries(WEAVER_POLYMERASES);
	polymeraseEntries.forEach(([key, polymerase]) => {
		const option = document.createElement("option");
		option.value = key;
		option.textContent = polymerase.name;
		selectElement.appendChild(option);
	});
	if (!polymeraseEntries.length) {
		const noPolymeraseOption = document.createElement("option");
		noPolymeraseOption.value = WEAVER_NO_POLYMERASE_OPTION;
		noPolymeraseOption.textContent = "No polymerase";
		selectElement.appendChild(noPolymeraseOption);
	}
	selectElement.value = WEAVER_POLYMERASES[weaverSelectedPolymeraseKey] ?
		weaverSelectedPolymeraseKey :
		firstPolymeraseKey() || WEAVER_NO_POLYMERASE_OPTION;
	weaverSelectedPolymeraseKey = selectElement.value === WEAVER_NO_POLYMERASE_OPTION ?
		"" :
		selectElement.value;
}

function expandAmpliconsForOve(annotations) {
	const sequenceLength = editorState.sequenceData && editorState.sequenceData.sequence ?
		editorState.sequenceData.sequence.length :
		0;
	if (!sequenceLength) {
		return annotations || [];
	}
	return (annotations || []).flatMap((annotation) => {
		const start = Number(annotation.start);
		const end = Number(annotation.end);
		const isCircular = annotation.overlapsSelf || end < start;
		if (!isCircular) {
			return [annotation];
		}
		return [
			{
				...annotation,
				id: annotation.id + "-segment-a",
				end: sequenceLength - 1,
				overlapsSelf: false,
				notes: {
					...(annotation.notes || {}),
					weaver_parent_amplicon_id: [annotation.id],
					weaver_circular_segment: ["1"]
				}
			},
			{
				...annotation,
				id: annotation.id + "-segment-b",
				start: 0,
				overlapsSelf: false,
				notes: {
					...(annotation.notes || {}),
					weaver_parent_amplicon_id: [annotation.id],
					weaver_circular_segment: ["2"]
				}
			}
		];
	});
}

function moduloForWeaver(value, length) {
	return ((Number(value) % length) + length) % length;
}

function rotateWeaverRange(range, rotateTo, sequenceLength) {
	if (!range || typeof range.start === "undefined" || typeof range.end === "undefined") {
		return range;
	}
	const rotatedRange = {
		...range,
		start: moduloForWeaver(Number(range.start) - rotateTo, sequenceLength),
		end: moduloForWeaver(Number(range.end) - rotateTo, sequenceLength)
	};
	if (Array.isArray(range.locations)) {
		rotatedRange.locations = range.locations.map((location) => rotateWeaverRange(location, rotateTo, sequenceLength));
	}
	return rotatedRange;
}

function rotateWeaverAnnotations(annotations, rotateTo, sequenceLength) {
	return (annotations || []).map((annotation) => rotateWeaverRange(annotation, rotateTo, sequenceLength));
}

function rotateWeaverSequenceData(sequenceData, rotateTo) {
	const sequence = sequenceData && sequenceData.sequence ? String(sequenceData.sequence) : "";
	const sequenceLength = sequence.length;
	if (!sequenceLength || !rotateTo) {
		return cloneForWeaver(sequenceData);
	}
	const rotatedSequenceData = cloneForWeaver(sequenceData);
	rotatedSequenceData.sequence = sequence.slice(rotateTo) + sequence.slice(0, rotateTo);
	[
		"features",
		"parts",
		"primers",
		"cutsites",
		"translations",
		"orfs",
		"warnings",
		"assemblyPieces"
	].forEach((key) => {
		if (Array.isArray(rotatedSequenceData[key])) {
			rotatedSequenceData[key] = rotateWeaverAnnotations(rotatedSequenceData[key], rotateTo, sequenceLength);
		}
	});
	return rotatedSequenceData;
}

function weaverSequenceDataForSave(sequenceDataToSave) {
	const originalSequenceData = weaverDigestOriginalSequenceData || weaverAmpliconOriginalSequenceData;
	return originalSequenceData ? cloneForWeaver(originalSequenceData) : sequenceDataToSave;
}

function restoreDigestOriginalSequenceData() {
	if (!weaverDigestOriginalSequenceData) {
		return;
	}
	editorState.sequenceData = cloneForWeaver(weaverDigestOriginalSequenceData);
	weaverDigestOriginalSequenceData = null;
	weaverDigestRotation = 0;
}

function restoreAmpliconOriginalSequenceData() {
	if (!weaverAmpliconOriginalSequenceData) {
		return;
	}
	editorState.sequenceData = cloneForWeaver(weaverAmpliconOriginalSequenceData);
	weaverAmpliconOriginalSequenceData = null;
	weaverAmpliconRotation = 0;
}

function restoreWeaverVisualOriginalSequenceData() {
	restoreDigestOriginalSequenceData();
	restoreAmpliconOriginalSequenceData();
}

function digestRotationForResult(result) {
	const sites = result && result.cut_sites ? result.cut_sites : [];
	return sites.length > 1 ? Number(sites[0].position) : 0;
}

function digestDisplaySites(result) {
	const sequenceLength = editorState.sequenceData && editorState.sequenceData.sequence ?
		editorState.sequenceData.sequence.length :
		0;
	const rotation = weaverDigestRotation || 0;
	return (result.cut_sites || []).map((site) => ({
		...site,
		position: sequenceLength ? moduloForWeaver(Number(site.position) - rotation, sequenceLength) : Number(site.position)
	})).sort((left, right) => left.position - right.position);
}

function originalDigestSelection(selection) {
	const originalSequenceData = weaverDigestOriginalSequenceData || weaverAmpliconOriginalSequenceData;
	const originalSequence = originalSequenceData && originalSequenceData.sequence ?
		String(originalSequenceData.sequence) :
		"";
	const sequenceLength = originalSequence.length || (
		editorState.sequenceData && editorState.sequenceData.sequence ? editorState.sequenceData.sequence.length : 0
	);
	const rotation = weaverDigestRotation || weaverAmpliconRotation || 0;
	if (!sequenceLength || !rotation) {
		return selection;
	}
	return {
		...selection,
		start: moduloForWeaver(Number(selection.start) + rotation, sequenceLength),
		end: moduloForWeaver(Number(selection.end) + rotation, sequenceLength),
		overlapsSelf: Boolean(selection.overlapsSelf || Number(selection.start) > Number(selection.end))
	};
}

function ampliconRotationForAnnotations(annotations) {
	if (!annotations || annotations.length !== 1) {
		return 0;
	}
	const annotation = annotations[0];
	const start = Number(annotation.start);
	const end = Number(annotation.end);
	return annotation.overlapsSelf || end < start ? start : 0;
}

function rotatedAmpliconAnnotations(annotations, rotation) {
	const sequenceLength = editorState.sequenceData && editorState.sequenceData.sequence ?
		editorState.sequenceData.sequence.length :
		0;
	if (!sequenceLength || !rotation) {
		return annotations || [];
	}
	return (annotations || []).map((annotation) => {
		const rotated = rotateWeaverRange(annotation, rotation, sequenceLength);
		rotated.overlapsSelf = Number(rotated.start) > Number(rotated.end);
		return rotated;
	});
}

function removeWeaverDigestAnnotations(sequenceData) {
	sequenceData.cutsites = (sequenceData.cutsites || []).filter((cutsite) => {
		return !(cutsite.id && String(cutsite.id).startsWith("weaver-digest-"));
	});
	sequenceData.parts = (sequenceData.parts || []).filter((part) => {
		return !(part.id && String(part.id).startsWith("weaver-digest-"));
	});
	return sequenceData;
}

function removeWeaverTemporaryFragments(sequenceData) {
	sequenceData = removeWeaverAmplicons(sequenceData);
	sequenceData = removeWeaverDigestAnnotations(sequenceData);
	return sequenceData;
}

function clearWeaverTemporarySelection() {
	editorState.selectionLayer = {
		start: -1,
		end: -1,
		forceUpdate: Math.random()
	};
	weaverLatestSelection = null;
}

function digestRegionRows() {
	return Array.from(document.querySelectorAll("#digest-regions .weaver-digest-region-row"));
}

function ampliconRegionRows() {
	return Array.from(document.querySelectorAll("#amplicon-regions .weaver-amplicon-region-row"));
}

function ampliconPrimerBindingRegionRows(direction) {
	return Array.from(document.querySelectorAll(
		"#amplicon-" + direction + "-binding-regions .weaver-amplicon-region-row"
	));
}

function digestRequiredEnzymeButtons() {
	return Array.from(document.querySelectorAll(".weaver-digest-enzyme-option"));
}

function selectedDigestRequiredEnzymes() {
	return digestRequiredEnzymeButtons()
		.filter((button) => button.classList.contains("is-selected"))
		.map((button) => button.dataset.enzyme)
		.filter(Boolean);
}

function setDigestRequiredEnzymeSelected(enzymeName, isSelected) {
	digestRequiredEnzymeButtons().forEach((button) => {
		if (button.dataset.enzyme !== enzymeName) {
			return;
		}
		button.classList.toggle("is-selected", Boolean(isSelected));
		button.setAttribute("aria-selected", isSelected ? "true" : "false");
	});
	renderDigestRequiredEnzymeChips();
}

function renderDigestRequiredEnzymeChips() {
	const chipBar = document.getElementById("digest-required-enzyme-chips");
	if (!chipBar) {
		return;
	}
	chipBar.innerHTML = "";
	const selectedButtons = digestRequiredEnzymeButtons().filter((button) => button.classList.contains("is-selected"));
	if (!selectedButtons.length) {
		const empty = document.createElement("span");
		empty.className = "weaver-digest-required-enzyme-empty";
		empty.textContent = "Any enzyme";
		chipBar.appendChild(empty);
		return;
	}
	selectedButtons.forEach((optionButton) => {
		const chip = document.createElement("span");
		chip.className = "weaver-digest-required-enzyme-chip";
		chip.textContent = optionButton.textContent;
		const removeButton = document.createElement("button");
		removeButton.type = "button";
		removeButton.title = "Remove " + optionButton.textContent;
		removeButton.setAttribute("aria-label", "Remove " + optionButton.textContent);
		removeButton.innerHTML = '<i class="bi bi-x"></i>';
		removeButton.onclick = () => setDigestRequiredEnzymeSelected(optionButton.dataset.enzyme, false);
		chip.appendChild(removeButton);
		chipBar.appendChild(chip);
	});
}

function setDigestRequiredEnzymePickerOpen(isOpen) {
	const picker = document.getElementById("digest-required-enzymes");
	const list = document.getElementById("digest-required-enzyme-list");
	const toggle = document.getElementById("digest-required-enzyme-toggle");
	if (!picker || !list || !toggle) {
		return;
	}
	picker.classList.toggle("is-open", Boolean(isOpen));
	list.classList.toggle("is-collapsed", !isOpen);
	toggle.setAttribute("aria-expanded", isOpen ? "true" : "false");
}

function mountDigestRequiredEnzymePicker() {
	const toggle = document.getElementById("digest-required-enzyme-toggle");
	const list = document.getElementById("digest-required-enzyme-list");
	if (toggle && list) {
		toggle.onclick = () => setDigestRequiredEnzymePickerOpen(list.classList.contains("is-collapsed"));
	}
	digestRequiredEnzymeButtons().forEach((button) => {
		button.onclick = () => {
			setDigestRequiredEnzymeSelected(
				button.dataset.enzyme,
				!button.classList.contains("is-selected")
			);
		};
	});
	renderDigestRequiredEnzymeChips();
	setDigestRequiredEnzymePickerOpen(false);
}

function currentOveSelectionForDigest() {
	const selection = weaverLatestSelection || (editorState && editorState.selectionLayer ? editorState.selectionLayer : null);
	if (!selection || Array.isArray(selection) || selection.start < 0 || selection.end < 0) {
		return null;
	}
	const originalSelection = originalDigestSelection(selection);
	return {
		start: Number(originalSelection.start) + 1,
		end: Number(originalSelection.end) + 1
	};
}

function currentOveSelectionForAmplicon() {
	const selection = weaverLatestSelection || (editorState && editorState.selectionLayer ? editorState.selectionLayer : null);
	if (!selection || Array.isArray(selection) || selection.start < 0 || selection.end < 0) {
		return null;
	}
	const originalSelection = originalDigestSelection(selection);
	return {
		start: Number(originalSelection.start) + 1,
		end: Number(originalSelection.end) + 1
	};
}

function digestRegionExists(startValue, endValue) {
	return digestRegionRows().some((row) => {
		const inputs = row.querySelectorAll("input");
		const rowStart = inputs[0] ? String(inputs[0].value || "").trim() : "";
		const rowEnd = inputs[1] ? String(inputs[1].value || "").trim() : "";
		return rowStart === String(startValue) && rowEnd === String(endValue);
	});
}

function ampliconRegionExists(startValue, endValue) {
	return ampliconRegionRows().some((row) => {
		const inputs = row.querySelectorAll("input");
		const rowStart = inputs[0] ? String(inputs[0].value || "").trim() : "";
		const rowEnd = inputs[1] ? String(inputs[1].value || "").trim() : "";
		return rowStart === String(startValue) && rowEnd === String(endValue);
	});
}

function ampliconPrimerBindingRegionExists(direction, startValue, endValue) {
	return ampliconPrimerBindingRegionRows(direction).some((row) => {
		const inputs = row.querySelectorAll("input");
		const rowStart = inputs[0] ? String(inputs[0].value || "").trim() : "";
		const rowEnd = inputs[1] ? String(inputs[1].value || "").trim() : "";
		return rowStart === String(startValue) && rowEnd === String(endValue);
	});
}

function addDigestRegion(startValue, endValue) {
	const container = document.getElementById("digest-regions");
	if (!container) {
		return;
	}
	const row = document.createElement("div");
	row.className = "weaver-digest-region-row";
	const startInput = document.createElement("input");
	startInput.className = "form-control form-control-sm";
	startInput.type = "text";
	startInput.inputMode = "numeric";
	startInput.pattern = "[0-9]*";
	startInput.placeholder = "Start";
	startInput.value = startValue || "";
	const endInput = document.createElement("input");
	endInput.className = "form-control form-control-sm";
	endInput.type = "text";
	endInput.inputMode = "numeric";
	endInput.pattern = "[0-9]*";
	endInput.placeholder = "End";
	endInput.value = endValue || "";
	const useSelectionButton = document.createElement("button");
	useSelectionButton.type = "button";
	useSelectionButton.className = "btn btn-sm btn-outline-primary";
	useSelectionButton.title = "Use active map selection";
	useSelectionButton.setAttribute("aria-label", "Use active map selection");
	useSelectionButton.innerHTML = '<i class="bi bi-bounding-box"></i>';
	useSelectionButton.onclick = () => {
		const selection = currentOveSelectionForDigest();
		if (!selection) {
			setDigestStatus("No active map selection is available.", true);
			return;
		}
		startInput.value = selection.start;
		endInput.value = selection.end;
	};
	const removeButton = document.createElement("button");
	removeButton.type = "button";
	removeButton.className = "btn btn-sm btn-outline-secondary";
	removeButton.title = "Remove region";
	removeButton.setAttribute("aria-label", "Remove region");
	removeButton.innerHTML = '<i class="bi bi-trash"></i>';
	removeButton.onclick = () => row.remove();
	row.appendChild(startInput);
	row.appendChild(endInput);
	row.appendChild(useSelectionButton);
	row.appendChild(removeButton);
	container.appendChild(row);
}

function addAmpliconRegion(startValue, endValue) {
	const container = document.getElementById("amplicon-regions");
	if (!container) {
		return;
	}
	const row = document.createElement("div");
	row.className = "weaver-digest-region-row weaver-amplicon-region-row";
	const startInput = document.createElement("input");
	startInput.className = "form-control form-control-sm";
	startInput.type = "text";
	startInput.inputMode = "numeric";
	startInput.pattern = "[0-9]*";
	startInput.placeholder = "Start";
	startInput.value = startValue || "";
	const endInput = document.createElement("input");
	endInput.className = "form-control form-control-sm";
	endInput.type = "text";
	endInput.inputMode = "numeric";
	endInput.pattern = "[0-9]*";
	endInput.placeholder = "End";
	endInput.value = endValue || "";
	const useSelectionButton = document.createElement("button");
	useSelectionButton.type = "button";
	useSelectionButton.className = "btn btn-sm btn-outline-primary";
	useSelectionButton.title = "Use active map selection";
	useSelectionButton.setAttribute("aria-label", "Use active map selection");
	useSelectionButton.innerHTML = '<i class="bi bi-bounding-box"></i>';
	useSelectionButton.onclick = () => {
		const selection = currentOveSelectionForAmplicon();
		if (!selection) {
			setPrimerMatchStatus("No active map selection is available.", true);
			return;
		}
		startInput.value = selection.start;
		endInput.value = selection.end;
	};
	const removeButton = document.createElement("button");
	removeButton.type = "button";
	removeButton.className = "btn btn-sm btn-outline-secondary";
	removeButton.title = "Remove region";
	removeButton.setAttribute("aria-label", "Remove region");
	removeButton.innerHTML = '<i class="bi bi-trash"></i>';
	removeButton.onclick = () => row.remove();
	row.appendChild(startInput);
	row.appendChild(endInput);
	row.appendChild(useSelectionButton);
	row.appendChild(removeButton);
	container.appendChild(row);
}

function addAmpliconPrimerBindingRegion(direction, startValue, endValue) {
	const container = document.getElementById("amplicon-" + direction + "-binding-regions");
	if (!container) {
		return;
	}
	const row = document.createElement("div");
	row.className = "weaver-digest-region-row weaver-amplicon-region-row";
	const startInput = document.createElement("input");
	startInput.className = "form-control form-control-sm";
	startInput.type = "text";
	startInput.inputMode = "numeric";
	startInput.pattern = "[0-9]*";
	startInput.placeholder = "Start";
	startInput.value = startValue || "";
	const endInput = document.createElement("input");
	endInput.className = "form-control form-control-sm";
	endInput.type = "text";
	endInput.inputMode = "numeric";
	endInput.pattern = "[0-9]*";
	endInput.placeholder = "End";
	endInput.value = endValue || "";
	const useSelectionButton = document.createElement("button");
	useSelectionButton.type = "button";
	useSelectionButton.className = "btn btn-sm btn-outline-primary";
	useSelectionButton.title = "Use active map selection";
	useSelectionButton.setAttribute("aria-label", "Use active map selection");
	useSelectionButton.innerHTML = '<i class="bi bi-bounding-box"></i>';
	useSelectionButton.onclick = () => {
		const selection = currentOveSelectionForAmplicon();
		if (!selection) {
			setPrimerMatchStatus("No active map selection is available.", true);
			return;
		}
		startInput.value = selection.start;
		endInput.value = selection.end;
	};
	const removeButton = document.createElement("button");
	removeButton.type = "button";
	removeButton.className = "btn btn-sm btn-outline-secondary";
	removeButton.title = "Remove region";
	removeButton.setAttribute("aria-label", "Remove region");
	removeButton.innerHTML = '<i class="bi bi-trash"></i>';
	removeButton.onclick = () => row.remove();
	row.appendChild(startInput);
	row.appendChild(endInput);
	row.appendChild(useSelectionButton);
	row.appendChild(removeButton);
	container.appendChild(row);
}

function digestRequestParams() {
	const params = new URLSearchParams({
		min_fragments: document.getElementById("digest-min-fragments").value || "2",
		max_fragments: document.getElementById("digest-max-fragments").value || "3",
		min_band_difference_bp: document.getElementById("digest-min-band-difference").value || "500",
		min_fragment_size_bp: document.getElementById("digest-min-fragment-size").value || "500",
		max_enzymes: document.getElementById("digest-max-enzymes").value || "2",
		limit: "10"
	});
	const requiredEnzymes = selectedDigestRequiredEnzymes();
	if (requiredEnzymes.length) {
		params.set("required_enzymes", requiredEnzymes.join(","));
	}
	const regions = digestRegionRows().map((row) => {
		const inputs = row.querySelectorAll("input");
		return {
			start: inputs[0] ? inputs[0].value : "",
			end: inputs[1] ? inputs[1].value : ""
		};
	}).filter((region) => String(region.start).trim() && String(region.end).trim());
	params.set("regions", JSON.stringify(regions));
	return params;
}

function ampliconRequestParams() {
	const params = new URLSearchParams({
		max_tm_diff: "5",
		non_overlapping: true
	});
	const minSize = document.getElementById("amplicon-min-size");
	const maxSize = document.getElementById("amplicon-max-size");
	const minSizeValue = minSize ? String(minSize.value || "").trim() : "";
	const maxSizeValue = maxSize ? String(maxSize.value || "").trim() : "";
	if (minSizeValue) {
		params.set("min_size", minSizeValue);
	}
	if (maxSizeValue) {
		params.set("max_size", maxSizeValue);
	}
	const primerIds = [
		document.getElementById("amplicon-primer-id-1"),
		document.getElementById("amplicon-primer-id-2")
	].map((input) => input ? String(input.value || "").trim() : "").filter(Boolean);
	if (primerIds.length) {
		params.set("primer_ids", primerIds.join(","));
	}
	const regions = ampliconRegionRows().map((row) => {
		const inputs = row.querySelectorAll("input");
		return {
			start: inputs[0] ? inputs[0].value : "",
			end: inputs[1] ? inputs[1].value : ""
		};
	}).filter((region) => String(region.start).trim() && String(region.end).trim());
	params.set("regions", JSON.stringify(regions));
	["fwd", "rev"].forEach((direction) => {
		const bindingRegions = ampliconPrimerBindingRegionRows(direction).map((row) => {
			const inputs = row.querySelectorAll("input");
			return {
				start: inputs[0] ? inputs[0].value : "",
				end: inputs[1] ? inputs[1].value : ""
			};
		}).filter((region) => String(region.start).trim() && String(region.end).trim());
		params.set(direction + "_binding_regions", JSON.stringify(bindingRegions));
	});
	return params;
}

function setDigestStatus(message, isError) {
	const summary = document.getElementById("ove-digest-panel-summary");
	if (!summary) {
		return;
	}
	summary.textContent = message || "";
	summary.classList.toggle("text-danger", Boolean(isError));
	summary.classList.toggle("text-secondary", !isError);
}

function formatDigestList(values) {
	return (values || []).join(", ");
}

function digestBufferLabel(result) {
	const buffer = result.best_buffer;
	if (!buffer) {
		return "buffer: none";
	}
	const activities = Object.entries(buffer.activities || {}).map(([enzyme, activity]) => enzyme + " " + activity + "%");
	return buffer.name + " (" + activities.join(", ") + ")";
}

function digestRegionLabel(result) {
	const regions = result.regions || [];
	if (!regions.length) {
		return "regions: any";
	}
	const covered = regions.filter((region) => region.covered).length;
	return "regions: " + covered + "/" + regions.length;
}

function digestViolationsText(result) {
	return (result.violations || []).map((violation) => violation.message).join("; ");
}

function digestCutsiteAnnotations(result) {
	return (result.cut_sites || []).map((site, index) => ({
		id: result.id + "-cut-" + index,
		name: site.enzymes.join("+") + " cut " + site.position_label,
		start: site.position,
		end: site.position,
		annotationTypePlural: "cutsites",
		color: "#d63384",
		notes: {
			weaver_digest: ["true"],
			enzymes: site.enzymes,
			position: [String(site.position_label)]
		}
	}));
}

function digestFragmentAnnotations(result) {
	const sequenceLength = editorState.sequenceData && editorState.sequenceData.sequence ?
		editorState.sequenceData.sequence.length :
		0;
	const sites = digestDisplaySites(result);
	if (!sequenceLength || sites.length < 2) {
		return [];
	}
	return sites.flatMap((site, index) => {
		const next = sites[(index + 1) % sites.length];
		const start = site.position;
		const end = next.position === 0 ? sequenceLength - 1 : next.position - 1;
		const fragmentSize = String(result.fragments_map_order[index] || "");
		const baseAnnotation = {
			name: "Digest fragment " + (result.fragments_map_order[index] || "") + " bp",
			annotationTypePlural: "parts",
			color: "#20c997",
			notes: {
				weaver_digest: ["true"],
				fragment_size: [fragmentSize]
			}
		};
		if (next.position > 0 && next.position <= site.position) {
			return [
				{
					...baseAnnotation,
					id: result.id + "-fragment-" + index + "-segment-a",
					start,
					end: sequenceLength - 1,
					overlapsSelf: false,
					notes: {
						...baseAnnotation.notes,
						weaver_circular_segment: ["1"]
					}
				},
				{
					...baseAnnotation,
					id: result.id + "-fragment-" + index + "-segment-b",
					start: 0,
					end,
					overlapsSelf: false,
					notes: {
						...baseAnnotation.notes,
						weaver_circular_segment: ["2"]
					}
				}
			];
		}
		return [{
			...baseAnnotation,
			id: result.id + "-fragment-" + index,
			start,
			end,
			overlapsSelf: false
		}];
	});
}

function showDigestResult(result) {
	if (!editor || !editorState.sequenceData || !result) {
		return;
	}
	restoreWeaverVisualOriginalSequenceData();
	editorState.sequenceData = removeWeaverTemporaryFragments(editorState.sequenceData);
	weaverVisibleAmpliconIds = new Set();
	const rotation = digestRotationForResult(result);
	if (rotation) {
		weaverDigestOriginalSequenceData = cloneForWeaver(editorState.sequenceData);
		weaverDigestRotation = rotation;
		editorState.sequenceData = rotateWeaverSequenceData(editorState.sequenceData, rotation);
	}
	editorState.sequenceData = removeWeaverTemporaryFragments(editorState.sequenceData);
	editorState.sequenceData.parts = [
		...(editorState.sequenceData.parts || []),
		...digestFragmentAnnotations(result)
	];
	weaverVisibleDigestId = result.id;
	editorState.annotationVisibility = {
		...(editorState.annotationVisibility || {}),
		parts: true
	};
	if ((result.cut_sites || []).length) {
		const firstSite = result.cut_sites[0];
		editorState.selectionLayer = {
			start: firstSite.position,
			end: firstSite.position,
			overlapsSelf: false,
			forceUpdate: Math.random(),
			color: "rgba(214, 51, 132, 0.28)"
		};
	}
	editor.updateEditor(editorState);
	setTimeout(mountWeaverPrimerToolbar, 0);
	renderDigestPanel(weaverDigestResults);
	setDigestStatus(result.enzymes.join(" + ") + " shown on map.");
}

function clearDigestAnnotations() {
	if (!editor || !editorState.sequenceData) {
		return;
	}
	restoreDigestOriginalSequenceData();
	editorState.sequenceData = removeWeaverDigestAnnotations(editorState.sequenceData);
	weaverVisibleDigestId = "";
	clearWeaverTemporarySelection();
	editor.updateEditor(editorState);
	setTimeout(mountWeaverPrimerToolbar, 0);
	renderDigestPanel(weaverDigestResults);
	setDigestStatus("Digest annotations cleared.");
}

function renderDigestPanel(results) {
	const list = document.getElementById("ove-digest-panel-list");
	if (!list) {
		return;
	}
	list.innerHTML = "";
	if (!results || !results.length) {
		const empty = document.createElement("div");
		empty.className = "small text-secondary py-2";
		empty.textContent = "No digest results yet.";
		list.appendChild(empty);
		return;
	}
	results.forEach((result) => {
		const item = document.createElement("div");
		item.className = "weaver-amplicon-item weaver-digest-item";
		item.classList.toggle("is-visible", weaverVisibleDigestId === result.id);
		const title = document.createElement("div");
		title.className = "weaver-amplicon-item-name";
		title.textContent = result.enzymes.join(" + ");
		item.appendChild(title);
		const meta = document.createElement("div");
		meta.className = "weaver-amplicon-item-meta";
		[
			result.status,
			digestBufferLabel(result),
			result.temperature ? "temp " + result.temperature + " °C" : "temp " + result.temperature_status,
			result.cut_count + " cuts",
			result.fragment_count + " fragments",
			"map " + formatDigestList(result.fragments_map_order) + " bp",
			"bands " + formatDigestList(result.fragments_by_size) + " bp",
			digestRegionLabel(result)
		].forEach((text) => {
			const tag = document.createElement("span");
			tag.textContent = text;
			meta.appendChild(tag);
		});
		const cutsChip = document.createElement("span");
		cutsChip.className = "weaver-digest-cuts-chip";
		cutsChip.textContent = "Cuts: " + (result.cut_sites || []).map((site) => site.position_label + " " + site.enzymes.join("+")).join("; ");
		meta.appendChild(cutsChip);
		item.appendChild(meta);
		if (!result.exact && (result.violations || []).length) {
			const violations = document.createElement("div");
			violations.className = "weaver-digest-violations";
			violations.textContent = digestViolationsText(result);
			item.appendChild(violations);
		}
		const actions = document.createElement("div");
		actions.className = "weaver-amplicon-item-actions";
		const button = document.createElement("button");
		button.type = "button";
		button.className = "btn btn-sm btn-outline-primary";
		button.innerHTML = '<i class="bi bi-eye"></i>';
		button.title = "Show this digest on the map";
		button.setAttribute("aria-label", "Show this digest on the map");
		button.onclick = () => showDigestResult(result);
		actions.appendChild(button);
		item.appendChild(actions);
		list.appendChild(item);
	});
}

async function loadRestrictionDigests() {
	if (typeof plasmid_restriction_digest_path === "undefined") {
		setDigestStatus("Restriction digest endpoint is not available.", true);
		return;
	}
	setDigestStatus("Searching restriction digests...");
	try {
		const response = await fetch(plasmid_restriction_digest_path + "?" + digestRequestParams().toString(), {
			headers: {
				"Accept": "application/json"
			}
		});
		const data = await response.json();
		if (!response.ok) {
			setDigestStatus(data.error || "Unable to search restriction digests.", true);
			return;
		}
		weaverDigestResults = data.results || [];
		weaverVisibleDigestId = "";
		renderDigestPanel(weaverDigestResults);
		setDigestStatus((data.exact_count || 0) + " exact, " + Math.max(0, (data.count || 0) - (data.exact_count || 0)) + " closest.");
	} catch (error) {
		setDigestStatus("Unable to search restriction digests.", true);
	}
}

function showAmpliconAnnotations(annotations, message) {
	if (!editor || !editorState.sequenceData) {
		return;
	}
	restoreWeaverVisualOriginalSequenceData();
	editorState.sequenceData = removeWeaverTemporaryFragments(editorState.sequenceData);
	weaverVisibleDigestId = "";
	const rotation = ampliconRotationForAnnotations(annotations || []);
	if (rotation) {
		weaverAmpliconOriginalSequenceData = cloneForWeaver(editorState.sequenceData);
		weaverAmpliconRotation = rotation;
		editorState.sequenceData = rotateWeaverSequenceData(editorState.sequenceData, rotation);
	}
	const displayAnnotations = rotatedAmpliconAnnotations(annotations || [], rotation);
	editorState.sequenceData = removeWeaverTemporaryFragments(editorState.sequenceData);
	editorState.sequenceData.parts = [
		...(editorState.sequenceData.parts || []),
		...expandAmpliconsForOve(displayAnnotations)
	];
	weaverVisibleAmpliconIds = new Set((annotations || []).map((annotation) => annotation.id));
	if (displayAnnotations.length === 1) {
		const selectionRange = rangeForAmpliconSelection(displayAnnotations[0]);
		const forceUpdate = Math.random();
		const selectionColor = "rgba(172, 104, 204, 0.28)";
		editorState.selectionLayer = Array.isArray(selectionRange) ?
			selectionRange.map((range) => ({
				...range,
				forceUpdate,
				color: selectionColor
			})) :
			{
				...selectionRange,
				forceUpdate,
				color: selectionColor
			};
		editorState.caretPosition = -1;
	}
	hideRestrictionEnzymeLayers();
	editor.updateEditor(editorState);
	setTimeout(mountWeaverPrimerToolbar, 0);
	renderAmpliconPanel(weaverAmpliconCandidates);
	renderDigestPanel(weaverDigestResults);
	if (message) {
		setPrimerMatchStatus(message);
	}
}

function clearAmpliconsFromMap() {
	if (!editor || !editorState.sequenceData) {
		return;
	}
	restoreAmpliconOriginalSequenceData();
	editorState.sequenceData = removeWeaverAmplicons(editorState.sequenceData);
	weaverVisibleAmpliconIds = new Set();
	clearWeaverTemporarySelection();
	editor.updateEditor(editorState);
	setTimeout(mountWeaverPrimerToolbar, 0);
	renderAmpliconPanel(weaverAmpliconCandidates);
	setPrimerMatchStatus("Amplicons cleared from map.");
}

function renderAmpliconPanel(candidates) {
	const list = document.getElementById("ove-amplicon-panel-list");
	const summary = document.getElementById("ove-amplicon-panel-summary");
	if (!list || !summary) {
		return;
	}
	list.innerHTML = "";
	summary.classList.remove("is-filter-fallback");
	summary.textContent = (candidates || []).length + " candidate primer pairs.";

	if (!candidates || candidates.length === 0) {
		const empty = document.createElement("div");
		empty.className = "small text-secondary py-2";
		empty.textContent = "No primer pairs match the current filters.";
		list.appendChild(empty);
		return;
	}

	candidates.forEach((amplicon) => {
		const item = document.createElement("div");
		item.className = "weaver-amplicon-item";
		item.classList.toggle("is-visible", weaverVisibleAmpliconIds.has(amplicon.id));

		const title = document.createElement("div");
		title.className = "weaver-amplicon-item-name";
		title.textContent = amplicon.name || "Amplicon";
			item.appendChild(title);

			const meta = document.createElement("div");
			meta.className = "weaver-amplicon-item-meta";
			[
				{text: valueFromNotes(amplicon, "product_size") + " bp"},
				{text: primerIdLabel(amplicon)},
				{text: ampliconCoordinateLabel(amplicon)},
				{text: annealingTemperatureLabel(amplicon), title: annealingTemperatureTooltip(amplicon), className: "weaver-amplicon-ta-chip"},
				{text: primerComplementarityLabel(amplicon), className: "weaver-amplicon-dimer-chip", title: "Show primer complementarity", onClick: () => togglePrimerComplementarityDetails(amplicon)},
				{text: "ΔTm " + valueFromNotes(amplicon, "tm_difference") + " °C"},
				{text: partialPrimerBindingLabel(amplicon), title: partialPrimerBindingTooltip(amplicon), className: "weaver-amplicon-partial-chip"},
				{text: extensionTimeLabel(amplicon)}
			].forEach((itemMeta) => {
				if (!itemMeta.text.trim()) {
					return;
				}
				const tag = document.createElement("span");
				tag.textContent = itemMeta.text;
				if (itemMeta.title) {
					tag.title = itemMeta.title;
				}
				if (itemMeta.className) {
					tag.className = itemMeta.className;
				}
				if (itemMeta.onClick) {
					tag.setAttribute("role", "button");
					tag.tabIndex = 0;
					tag.onclick = itemMeta.onClick;
					tag.onkeydown = (event) => {
						if (event.key === "Enter" || event.key === " ") {
							event.preventDefault();
							itemMeta.onClick();
						}
					};
				}
				meta.appendChild(tag);
			});
			item.appendChild(meta);

			const dimerExpanded = weaverExpandedDimerIds.has(amplicon.id);
			const warnings = primerComplementarityWarnings(amplicon);
			const alignment = primerComplementarityAlignment(amplicon);
			const primer3Structure = primer3DimerStructure(amplicon);
			if (dimerExpanded) {
				const detailBlock = document.createElement("div");
				detailBlock.className = "weaver-amplicon-dimer-detail";

				const summaryBlock = document.createElement("div");
				summaryBlock.className = "weaver-amplicon-item-warnings";
				summaryBlock.textContent = primerComplementaritySummary(amplicon);
				detailBlock.appendChild(summaryBlock);

				if (warnings.length) {
					const warningBlock = document.createElement("div");
					warningBlock.className = "weaver-amplicon-item-warnings";
					warningBlock.textContent = warnings.join("; ");
					detailBlock.appendChild(warningBlock);
				}
				if (alignment.length) {
					const alignmentLabel = document.createElement("div");
					alignmentLabel.className = "weaver-amplicon-dimer-label";
					alignmentLabel.textContent = "Weaver complementarity scan";
					detailBlock.appendChild(alignmentLabel);
					const alignmentBlock = document.createElement("pre");
					alignmentBlock.className = "weaver-amplicon-dimer-preview";
					alignmentBlock.textContent = alignment.join("\n");
					detailBlock.appendChild(alignmentBlock);
				}
				if (primer3Structure) {
					const structureLabel = document.createElement("div");
					structureLabel.className = "weaver-amplicon-dimer-label";
					structureLabel.textContent = "Primer3 ASCII structure";
					detailBlock.appendChild(structureLabel);
					const structureBlock = document.createElement("pre");
					structureBlock.className = "weaver-amplicon-dimer-preview";
					structureBlock.textContent = primer3Structure;
					detailBlock.appendChild(structureBlock);
				}
				item.appendChild(detailBlock);
			}

			const actions = document.createElement("div");
			actions.className = "weaver-amplicon-item-actions";
			const button = document.createElement("button");
			button.type = "button";
		button.className = "btn btn-sm btn-outline-primary";
			button.innerHTML = '<i class="bi bi-eye"></i>';
			button.title = "Show this amplicon on the map";
			button.setAttribute("aria-label", "Show this amplicon on the map");
			button.onclick = () => showAmpliconAnnotations([amplicon], "1 amplicon shown on map.");
			actions.appendChild(button);
			const copyButton = document.createElement("button");
			copyButton.type = "button";
			copyButton.className = "btn btn-sm btn-outline-secondary ms-1";
			copyButton.innerHTML = '<i class="bi bi-clipboard"></i>';
			copyButton.title = "Copy amplicon sequence";
			copyButton.setAttribute("aria-label", "Copy amplicon sequence");
			copyButton.onclick = () => copyAmpliconSequence(amplicon);
			actions.appendChild(copyButton);
			const blastButton = document.createElement("button");
			blastButton.type = "button";
			blastButton.className = "btn btn-sm btn-outline-success ms-1";
			blastButton.innerHTML = '<i class="bi bi-search"></i>';
			blastButton.title = "BLAST amplicon against local plasmids";
			blastButton.setAttribute("aria-label", "BLAST amplicon against local plasmids");
			blastButton.onclick = () => blastAmpliconSequence(amplicon);
			actions.appendChild(blastButton);
			item.appendChild(actions);

		list.appendChild(item);
	});
}

function makePrimerToolButton(className, title, html, onClick) {
	const button = document.createElement("button");
	button.type = "button";
	button.className = className;
	button.title = title;
	button.setAttribute("aria-label", title);
	button.innerHTML = html;
	button.onclick = onClick;
	return button;
}

function mountWeaverPrimerToolbar() {
	const existing = document.getElementById("ove-primer-tools");
	if (existing) {
		existing.remove();
	}
	const panelTab = document.getElementById("ove-amplicon-panel-tab");
	if (panelTab) {
		panelTab.onclick = () => {
			setAmpliconPanelOpen(true);
			if (!weaverAmpliconCandidatesLoaded) {
				loadPlasmidAmplicons();
			}
		};
	}
	const panelClose = document.getElementById("ove-amplicon-panel-close");
	if (panelClose) {
		panelClose.onclick = () => setAmpliconPanelOpen(false);
	}
	const panelClearMap = document.getElementById("ove-amplicon-panel-clear-map");
	if (panelClearMap) {
		panelClearMap.onclick = clearAmpliconsFromMap;
	}
	const ampliconParametersToggle = document.getElementById("amplicon-parameters-toggle");
	const ampliconParametersBody = document.getElementById("amplicon-parameters-body");
	const ampliconPanel = document.getElementById("ove-amplicon-panel");
	if (ampliconParametersToggle && ampliconParametersBody && ampliconPanel) {
		ampliconParametersToggle.onclick = () => {
			const isCollapsed = !ampliconParametersBody.classList.contains("is-collapsed");
			ampliconParametersBody.classList.toggle("is-collapsed", isCollapsed);
			ampliconPanel.classList.toggle("has-collapsed-parameters", isCollapsed);
			ampliconParametersToggle.setAttribute("aria-expanded", isCollapsed ? "false" : "true");
		};
	}
	const ampliconAddRegion = document.getElementById("amplicon-add-region");
	if (ampliconAddRegion) {
		ampliconAddRegion.onclick = () => {
			const selection = currentOveSelectionForAmplicon();
			if (!selection || ampliconRegionExists(selection.start, selection.end)) {
				addAmpliconRegion();
				return;
			}
			addAmpliconRegion(selection.start, selection.end);
		};
	}
	["fwd", "rev"].forEach((direction) => {
		const addBindingRegionButton = document.getElementById("amplicon-add-" + direction + "-binding-region");
		if (!addBindingRegionButton) {
			return;
		}
		addBindingRegionButton.onclick = () => {
			const selection = currentOveSelectionForAmplicon();
			if (!selection || ampliconPrimerBindingRegionExists(direction, selection.start, selection.end)) {
				addAmpliconPrimerBindingRegion(direction);
				return;
			}
			addAmpliconPrimerBindingRegion(direction, selection.start, selection.end);
		};
	});
	const ampliconRunSearch = document.getElementById("amplicon-run-search");
	if (ampliconRunSearch) {
		ampliconRunSearch.onclick = loadPlasmidAmplicons;
	}
	const digestPanelTab = document.getElementById("ove-digest-panel-tab");
	if (digestPanelTab) {
		digestPanelTab.onclick = () => setDigestPanelOpen(true);
	}
	const digestPanelClose = document.getElementById("ove-digest-panel-close");
	if (digestPanelClose) {
		digestPanelClose.onclick = () => setDigestPanelOpen(false);
	}
	const digestParametersToggle = document.getElementById("digest-parameters-toggle");
	const digestParametersBody = document.getElementById("digest-parameters-body");
	const digestPanel = document.getElementById("ove-digest-panel");
	if (digestParametersToggle && digestParametersBody && digestPanel) {
		digestParametersToggle.onclick = () => {
			const isCollapsed = !digestParametersBody.classList.contains("is-collapsed");
			digestParametersBody.classList.toggle("is-collapsed", isCollapsed);
			digestPanel.classList.toggle("has-collapsed-parameters", isCollapsed);
			digestParametersToggle.setAttribute("aria-expanded", isCollapsed ? "false" : "true");
		};
	}
	const digestAddRegion = document.getElementById("digest-add-region");
	if (digestAddRegion) {
		digestAddRegion.onclick = () => {
			const selection = currentOveSelectionForDigest();
			if (!selection || digestRegionExists(selection.start, selection.end)) {
				addDigestRegion();
				return;
			}
			addDigestRegion(
				selection.start,
				selection.end
			);
		};
	}
	const digestRunSearch = document.getElementById("digest-run-search");
	if (digestRunSearch) {
		digestRunSearch.onclick = loadRestrictionDigests;
	}
	const digestClearMap = document.getElementById("digest-clear-map");
	if (digestClearMap) {
		digestClearMap.onclick = clearDigestAnnotations;
	}
	mountDigestRequiredEnzymePicker();
	const polymeraseSelect = document.getElementById("ove-polymerase-select");
	if (polymeraseSelect) {
		populatePolymeraseSelect(polymeraseSelect);
		polymeraseSelect.onchange = () => {
			if (polymeraseSelect.value === WEAVER_NO_POLYMERASE_OPTION) {
				weaverSelectedPolymeraseKey = "";
				renderAmpliconPanel(weaverAmpliconCandidates);
				return;
			}
			weaverSelectedPolymeraseKey = polymeraseSelect.value;
			if (weaverSelectedPolymeraseKey === WEAVER_CUSTOM_POLYMERASE_KEY) {
				promptForPolymeraseSecondsPerKb(weaverSelectedPolymeraseKey);
				populatePolymeraseSelect(polymeraseSelect);
			}
			renderAmpliconPanel(weaverAmpliconCandidates);
		};
	}
}

function removeWeaverPrimerMatches(sequenceData) {
	const primers = sequenceData.primers || [];
	sequenceData.primers = primers.filter((primer) => {
		return !(primer.id && String(primer.id).startsWith("weaver-primer-"));
	});
	return sequenceData;
}

function removeWeaverAmplicons(sequenceData) {
	const parts = sequenceData.parts || [];
	sequenceData.parts = parts.filter((part) => {
		return !(part.id && String(part.id).startsWith("weaver-amplicon-"));
	});
	return sequenceData;
}

function clearPlasmidPrimerMatches() {
	if (!editor || !editorState.sequenceData) {
		return;
	}
	restoreWeaverVisualOriginalSequenceData();
	editorState.sequenceData = removeWeaverPrimerMatches(editorState.sequenceData);
	editorState.sequenceData = removeWeaverAmplicons(editorState.sequenceData);
	weaverAmpliconCandidates = [];
		weaverNonOverlappingAmplicons = [];
		weaverAmpliconCandidatesLoaded = false;
		weaverVisibleAmpliconIds = new Set();
		weaverExpandedDimerIds = new Set();
	editor.updateEditor(editorState);
	setTimeout(mountWeaverPrimerToolbar, 0);
	renderAmpliconPanel([]);
	setPrimerMatchStatus("Primer and amplicon matches cleared.");
}

function hideRestrictionEnzymeLayers() {
	editorState.annotationVisibility = {
		...(editorState.annotationVisibility || {}),
		primers: true,
		parts: true,
		cutsites: false,
		cutsiteLabels: false
	};
}

function isWeaverFocusedFeature(feature) {
	const type = String(feature.type || "").toLowerCase();
	const name = String(feature.name || feature.label || "").toLowerCase();
	const notes = feature.notes ? JSON.stringify(feature.notes).toLowerCase() : "";
	const searchableText = type + " " + name + " " + notes;
	const alwaysShowTypes = ["cds", "promoter", "terminator", "misc_recomb"];
	const recombinationKeywords = [
		"recomb",
		"scar",
		"conl",
		"conr",
		"cons",
		"lox",
		"frt",
		"attb",
		"attp",
		"attr",
		"attl"
	];
	return alwaysShowTypes.includes(type) || recombinationKeywords.some((keyword) => searchableText.includes(keyword));
}

function applyCdsPrimerView() {
	if (!editor || !editorState.sequenceData) {
		return;
	}
	const features = editorState.sequenceData.features || [];
	if (!editorState.sequenceData._weaverOriginalFeatures) {
		editorState.sequenceData._weaverOriginalFeatures = cloneForWeaver(features);
	}
	if (!editorState._weaverOriginalAnnotationVisibility) {
		editorState._weaverOriginalAnnotationVisibility = cloneForWeaver(editorState.annotationVisibility);
	}
	editorState.sequenceData.features = features.filter(isWeaverFocusedFeature);
	delete editorState.sequenceData.filteredFeatures;
	editorState.annotationVisibility = {
		...(editorState.annotationVisibility || {}),
		features: true,
		primers: true,
		parts: false,
		cutsites: false,
		cutsiteLabels: false,
		orfs: false,
		translations: false,
		orfTranslations: false,
		cdsFeatureTranslations: false,
		warnings: false,
		chromatogram: false,
		lineageLines: false,
		dnaColors: false,
		axis: true,
		axisNumbers: true,
		sequence: true,
		reverseSequence: true
	};
	editor.updateEditor(editorState);
	setTimeout(mountWeaverPrimerToolbar, 0);
	setPrimerMatchStatus(editorState.sequenceData.features.length + " functional features shown with primer matches.");
}

function restoreOveAnnotationView() {
	if (!editor || !editorState.sequenceData) {
		return;
	}
	restoreWeaverVisualOriginalSequenceData();
	if (editorState.sequenceData._weaverOriginalFeatures) {
		editorState.sequenceData.features = cloneForWeaver(editorState.sequenceData._weaverOriginalFeatures);
		delete editorState.sequenceData._weaverOriginalFeatures;
	}
	delete editorState.sequenceData.filteredFeatures;
	editorState.sequenceData = removeWeaverPrimerMatches(editorState.sequenceData);
	editorState.sequenceData = removeWeaverAmplicons(editorState.sequenceData);
	editorState.annotationVisibility = editorState._weaverDefaultAnnotationVisibility ?
		cloneForWeaver(editorState._weaverDefaultAnnotationVisibility) :
		cloneForWeaver(editorState._weaverOriginalAnnotationVisibility || {features: true});
	delete editorState._weaverOriginalAnnotationVisibility;
	editor.updateEditor(editorState);
	setTimeout(mountWeaverPrimerToolbar, 0);
	setPrimerMatchStatus("Annotation view restored.");
}

async function loadPlasmidAmplicons() {
	if (typeof plasmid_amplicon_matches_path === "undefined") {
		setPrimerMatchStatus("Amplicon endpoint is not available.", true);
		return;
	}
	if (!editor || !editorState.sequenceData) {
		setPrimerMatchStatus("OVE is not ready yet.", true);
		return;
	}
	setPrimerMatchStatus("Loading amplicons...");
	try {
		restoreWeaverVisualOriginalSequenceData();
		const params = ampliconRequestParams();
		const response = await fetch(plasmid_amplicon_matches_path + "?" + params.toString(), {
			headers: {
				"Accept": "application/json"
			}
		});
		const data = await response.json();
		if (!response.ok) {
			setPrimerMatchStatus(data.error || "Unable to load amplicons.", true);
			return;
		}
			weaverAmpliconCandidates = data.candidates || data.amplicons || [];
			weaverNonOverlappingAmplicons = data.amplicons || [];
			weaverAmpliconCandidatesLoaded = true;
			weaverExpandedDimerIds = new Set();
		if (editorState.sequenceData) {
			editorState.sequenceData = removeWeaverAmplicons(editorState.sequenceData);
			weaverVisibleAmpliconIds = new Set();
			editor.updateEditor(editorState);
			setTimeout(mountWeaverPrimerToolbar, 0);
		}
		renderAmpliconPanel(weaverAmpliconCandidates);
		const fallbackMessage = data.primer_filter_fallback ?
			" No exact pair found; showing combinations for each selected primer." : "";
		const statusMessage =
			(data.candidate_count || data.count || 0) + " candidate primer pairs loaded." + fallbackMessage;
		const panelSummary = document.getElementById("ove-amplicon-panel-summary");
		if (panelSummary) {
			panelSummary.classList.toggle("is-filter-fallback", Boolean(data.primer_filter_fallback));
			if (data.primer_filter_fallback) {
				panelSummary.textContent = statusMessage;
			}
		}
		setPrimerMatchStatus(statusMessage);
		setAmpliconPanelOpen(true);
	} catch (error) {
		setPrimerMatchStatus("Unable to load amplicons.", true);
	}
}

async function loadPlasmidPrimerMatches(options = {}) {
	if (typeof plasmid_primer_matches_path === "undefined") {
		setPrimerMatchStatus("Primer match endpoint is not available.", true);
		return;
	}
	if (!editor || !editorState.sequenceData) {
		setPrimerMatchStatus("OVE is not ready yet.", true);
		return;
	}
	setPrimerMatchStatus("Loading primer matches...");
	try {
		restoreWeaverVisualOriginalSequenceData();
		const response = await fetch(plasmid_primer_matches_path, {
			headers: {
				"Accept": "application/json"
			}
		});
		const data = await response.json();
		if (!response.ok) {
			setPrimerMatchStatus(data.error || "Unable to load primer matches.", true);
			return;
		}
		editorState.sequenceData = removeWeaverPrimerMatches(editorState.sequenceData);
		editorState.sequenceData.primers = [
			...(editorState.sequenceData.primers || []),
			...(data.primers || [])
		];
		hideRestrictionEnzymeLayers();
		editor.updateEditor(editorState);
		if (options.cdsOnly) {
			applyCdsPrimerView();
		} else {
			setTimeout(mountWeaverPrimerToolbar, 0);
			setPrimerMatchStatus((data.count || 0) + " primer matches added.");
		}
	} catch (error) {
		setPrimerMatchStatus("Unable to load primer matches.", true);
	}
}

async function loadCdsPrimerView() {
	await loadPlasmidPrimerMatches({cdsOnly: true});
}

async function sequenceToJson() {
    if(typeof sequence_file_contents == 'undefined') {return;}
    var jsonOutput = await window.bioParsers.anyToJson(sequence_file_contents);
    if(jsonOutput[0]['success']){
        editorState['sequenceData'] = jsonOutput[0]['parsedSequence'];
		editorState._weaverDefaultAnnotationVisibility = cloneForWeaver(editorState.annotationVisibility);
        editor = window.createVectorEditor(document.getElementById("ove-viewer"), editorProps);
        editor.updateEditor(editorState);
		stopCircularMapWheelPageScroll();
		setTimeout(mountWeaverPrimerToolbar, 0);
    } else {
        window.toastr.success("Error parsing plasmid sequence file.");
    }
}
loadCustomPolymeraseTime();
sequenceToJson();
