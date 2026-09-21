import React from 'react';
import './assembly_wizard.css';
import AssemblyStandards from './components/AssemblyStandards';
import * as ReactDOMServer from 'react-dom/server';
import $ from 'jquery';

const defaultState = {
    enzyme: null,
    receiver: null,
    inserts: [],
    complete: false,
    name: '',
    loading: false,
    parts: null,
    assembly_standard: 'ytk',
    apiError: '',
    receiverFilter: '',
    insertFilter: '',
    csrf_token: ''
}

function OHRender(ohSeq, assembly_standard) {
    let result = "Custom: " + ohSeq
    const standardOHS = AssemblyStandards[assembly_standard].ohs
    for (var key in standardOHS) {
        if (standardOHS[key].oh === ohSeq) {
            result = standardOHS[key].name + ": " + ohSeq
        }
    }
    return result
}

class PartRenderAssembly extends React.Component {
    render() {
        const part = this.props.part

        return <div className="alert border border-secondary">
            <p className="mb-0"><strong>{part.n}</strong> <a href={"/inventory/plasmid/" + part.i} target="_blank"
                                                             rel="noreferrer"><i
                className="bi bi-box-arrow-up-right"></i></a></p>
            <p className="mb-0"><span
                className="text-muted">{OHRender(part.o5, this.props.assembly_standard)} / {OHRender(part.o3, this.props.assembly_standard)}</span>
            </p>
        </div>
    }
}

class PartRender extends React.Component {
    render() {
        const part = this.props.part

        let className = ""
        if (this.props.active) {
            className = "table-primary"
        }

        let type_text = "Insert"
        if (part.t)
            type_text = "Receiver"

        return <tr className={className}>
            <td><span>{part.n}</span><a className="ms-2" href={"/inventory/plasmid/" + part.i} target="_blank"
                                        rel="noreferrer"><i
                class="bi bi-box-arrow-up-right"></i></a>
                <br/>
                <span className="small">OH 5 / </span><span
                    className="small text-secondary">{OHRender(part.o5, this.props.assembly_standard)}</span>
                <span className="small ms-2">OH 3 / </span><span
                    className="small text-secondary">{OHRender(part.o3, this.props.assembly_standard)}</span>
            </td>
            <td>
                <span>L{part.l}</span><span className="text-secondary"> - {type_text}</span>
            </td>
            <td className="text-break">
                {part.d}
            </td>
            <td>
                <button
                    type="button"
                    className="btn btn-outline-secondary"
                    value={part.i}
                    onClick={this.props.partHandler}>
                    {this.props.button_text}
                </button>
            </td>
        </tr>
    }
}

class EnzymeRender extends React.Component {
    render() {
        let button_class = "btn btn-outline-secondary"
        if (this.props.active) button_class = "btn btn-primary"

        return <button type="button" className={button_class + " me-2 mb-2"} name={"enzyme-" + this.props.value}
                       value={this.props.value} onClick={this.props.setEnzymeHandler}>{this.props.name}</button>
    }
}

class TableFilter extends React.Component {
    render() {
        return <div className="col-md-4">
            <div className="input-group">
                <input value={this.props.text} className="form-control" placeholder="Filter ..."
                       onChange={this.props.filterHandler}/>
                <button onClick={this.props.filterClearHandler} type="button" className="btn bg-transparent">
                    <i className="bi bi-x-circle"></i>
                </button>
            </div>
            <p></p>
        </div>
    }
}

class AssemblyResult extends React.Component {
    render() {
        const config = this.props.config
        const receiver = config.receiver

        let receiver_output = <div className="alert alert-info">No receiver selected</div>
        let complete_assembly = <div className="alert alert-info">Add parts to complete the assembly</div>
        let inserts_output = <div className="alert alert-info">No inserts selected</div>

        if (receiver) {
            receiver_output = []
            receiver_output.push(<p className="alert alert-primary">Backbone</p>)
            receiver_output.push(<PartRenderAssembly assembly_standard={config.assembly_standard} part={receiver}/>)
        }
        if (config.inserts.length) {
            if (config.complete) {
                complete_assembly = <div className="alert alert-success">
                    <div className="row">
                        <div className="col-12">
                            <span className={"d-inline-block"}>Assembly complete</span>
                            <button type="submit" form="form-wizard" value="Submit"
                                    className={"float-end md-2 btn btn-primary " + this.props.next_disabled}>{this.props.next_text}<i
            className="bi bi-arrow-right-circle ms-2"></i>  </button>
                        </div>
                    </div>
                </div>
            }
            inserts_output = []
            inserts_output.push(
                <div className="alert alert-primary">Inserts</div>
            )
            config.inserts.forEach((insert) => {
                inserts_output.push(<PartRenderAssembly assembly_standard={config.assembly_standard}
                                                        part={insert}/>)
            })
            inserts_output.push(
                <p>
                    <button type="button" className="btn btn-outline-danger"
                            onClick={this.props.removeLastPartHandler}>Remove last insert
                    </button>
                </p>
            )
        }

        return <div>
            <div>
                {complete_assembly}
            </div>
            <div>
                {receiver_output}
            </div>
            <div>
                {inserts_output}
            </div>
        </div>
    }
}

class PlasmidMap extends React.Component {
    constructor(props) {
        super(props)
        this.state = {
            canvas_size: 800
        }
        this.theIframe = React.createRef()
    }

    componentDidMount() {
        this.setState({
            canvas_size: this.theIframe.current.offsetWidth
        })
    }

    render() {
        const config = this.props.config

        if (config.receiver) {
            let length = config.receiver.len
            let insert_trackmarkers = []

            const colors = [
                'rgba(170,0,85,0.9)',
                'rgba(85,0,170,0.9)',
                'rgba(0,85,170,0.9)',
                'rgba(85,170,0,0.9)',
                'rgba(170,85,0,0.9)',
                'rgba(128,64,64,0.8)',
            ]

            let colori = 0

            if (config.inserts.length) {
                let currPos = 0
                config.inserts.forEach((insert) => {
                    length += insert.len
                    insert_trackmarkers.push(
                        <trackmarker start={currPos} end={currPos + insert.len} markerstyle={"fill:" + colors[colori]}>
                            <markerlabel vadjust="50" type="path"
                                         labelstyle={"font-size:12px;font-weight:400;fill:" + colors[colori]}
                                         text={insert.n}></markerlabel>
                        </trackmarker>
                    )
                    colori += 1
                    if (colori >= colors.length) colori = 0
                    currPos += insert.len
                })
            }

            let plasmid_name = this.props.name
            if (!plasmid_name) plasmid_name = "Set plasmid name"

            let canvas_size = this.state.canvas_size
            if (canvas_size > 600) canvas_size = 600

            let html = ReactDOMServer.renderToString(
                <html>
                <head>
                    <title>Plasmid Map</title>
                    <script src='/static/js/angularplasmid.complete.min.js'></script>
                </head>
                <body style={{'width': '100%', 'overflow': 'hidden', 'textAlign': 'center'}}>
                <plasmid id='p1' sequencelength={length} plasmidwidth={canvas_size + 50}
                         plasmidheight={canvas_size + 100}>
                    <plasmidtrack width="5" trackstyle='fill:#ccc;' radius={canvas_size / 2 - 40}></plasmidtrack>
                    <plasmidtrack trackstyle='fill:rgba(225,225,225,0.5);' radius={canvas_size / 2 - 50}>
                        <tracklabel labelstyle="font-size:20px;font-weight:600;fill:#666;"
                                    text={plasmid_name}></tracklabel>
                        <tracklabel labelstyle="font-size:12px;font-weight:200;fill:#999;" text={length + " bp"}
                                    vadjust="18"></tracklabel>
                        <trackscale className='smajor' style={{'strokeWidth': 1, stroke: "#888"}} ticksize='4'
                                    interval={200} labelstyle="font-size:9px;fill:#666" labelvadjust="15" showlabels='1'
                                    labelclass='sml'></trackscale>
                        <trackmarker start={length - config.receiver.len} end={length} markerstyle="fill:transparent">
                            <markerlabel vadjust="50" type="path" labelstyle="font-size:12px;fill:#888;font-weight:400"
                                         text={"Backbone: " + config.receiver.name}></markerlabel>
                        </trackmarker>
                        {insert_trackmarkers}
                    </plasmidtrack>
                </plasmid>
                </body>
                </html>
            )
            const iframe = ReactDOMServer.renderToString(
                <iframe title="Plasmid Map" width="100%" height={canvas_size + 100} srcDoc={html}></iframe>
            )
            return <div ref={this.theIframe} dangerouslySetInnerHTML={{__html: iframe}}></div>
        }

        return <div ref={this.theIframe} className="alert alert-info">Set a backbone to generate the plasmid map</div>
    }
}

class Pagination extends React.Component {
    render() {
        let next_button = <button className="btn btn-success float-end" onClick={() => {
            this.props.stepHandler(this.props.step + 1)
        }}>Next</button>
        if (this.props.step === 3)
            next_button = <button type="submit" form="form-wizard" value="Submit"
                                  className={"float-end me-2 btn btn-primary " + this.props.next_disabled}>{this.props.next_text}<i
                className="bi bi-arrow-right-circle ms-2"></i></button>
        let prev_button = ""
        if (this.props.step > 1)
            prev_button = <button className="btn btn-secondary float-end me-2" onClick={() => {
                this.props.stepHandler(this.props.step - 1)
            }}>Prev</button>
        let pre_hr = ""
        let post_hr = ""
        if (this.props.position === "t")
            post_hr = <hr/>
        if (this.props.position === "b")
            pre_hr = <hr/>
        return <div className="col-12">
            {pre_hr}
            <div className="flow_root">
                <a href={this.props.create_url} className="btn btn-outline-secondary float-start me-2">Skip wizard</a>
                {next_button}
                {prev_button}
            </div>
            {post_hr}
        </div>
    }
}

class App extends React.Component {
    constructor(props) {
        super(props)
        this.state = defaultState
    }

    addInsertHandler = (event) => {
        this.setState({
            insertFilter: "",
        })
        this.addPartHandler(event)
    }
    filterReceiverHandler = (event) => {
        this.setState({
            receiverFilter: event.target.value,
        })
    }
    filteriInsertHandler = (event) => {
        this.setState({
            insertFilter: event.target.value,
        })
    }
    filterReceiverClearHandler = (event) => {
        this.setState({
            receiverFilter: "",
        })
    }
    filteriInsertClearHandler = (event) => {
        this.setState({
            insertFilter: "",
        })
    }

    setStandardHandler = (event) => {
        this.setState({
            receiver: null,
            inserts: [],
            parts: [],
            complete: false,
            apiError: '',
            assembly_standard: event.target.value,
            enzyme: ''
        })
    }
    setName = (event) => {
        this.setState({
            name: event.target.value
        })
    }
    setDescription = (event) => {
        this.setState({
            description: event.target.value
        })
    }
    removeLastPartHandler = (event) => {
        let inserts = this.state.inserts
        inserts.pop()

        this.setState({
            inserts: inserts,
            complete: false
        })
    }
    addPartHandler = (event) => {
        let inserts = this.state.inserts
        let newPart = null
        this.state.parts.forEach((part) => {
            if (part.i === event.target.value) {
                newPart = part
            }
        })

        let complete = false
        if (newPart) {
            inserts.push(newPart)
            if(this.state.receiver){
                if (this.state.receiver.o5 === newPart.o3)
                    complete = true
            } else {
                console.log(inserts[0]);
                console.log(newPart);
                if (inserts[0].o5 === newPart.o3)
                    complete = true
            }
        }

        this.setState({
            inserts: inserts,
            complete: complete
        })
    }
    setReceiverHandler = (event) => {
        let receiver = null
        this.state.parts.forEach((part) => {
            if (part.i === event.target.value) {
                receiver = part
            }
        })
        let inserts = this.state.inserts
        if (this.state.receiver)
            if (receiver.oh5 !== this.state.receiver.oh5
                || receiver.oh3 !== this.state.receiver.oh3) {
                inserts = []
            }
        this.setState({
            receiver: receiver,
            inserts: inserts
        })
    }
    setEnzymeHandler = (event) => {
        if (event.target.value !== this.state.enzyme) {
            this.setState({
                receiver: null,
                inserts: [],
                parts: [],
                complete: false,
                loading: true,
                apiError: '',
                enzyme: event.target.value
            }, () => {
                this.apiCall()
            })
        }
    }

    apiCall() {
        const axios = require('axios');
        // const url = 'http://192.168.1.64:8000/inventory/api/parts/' + this.state.enzyme + '/' + this.state.assembly_standard + '/'
        const url = '/inventory/api/parts/' + this.state.enzyme + '/' + this.state.assembly_standard + '/'
        axios.get(url)
            .then((response) => {
                this.setState({
                    apiError: response.data.error,
                    parts: response.data.parts,
                    csrf_token: response.data.csrf_token,
                    loading: false
                })
            })
    }

    stepHandler(step) {
        $('#steps').children().removeClass('show')
        $('#steps > *:nth-child(' + (step) + ')').addClass('show')
        $('#pagination').children().removeClass('active')
        $('#pagination > *:nth-child(' + (step) + ')').addClass('active')
        $('#pagination a span:first-child').removeClass('bg-light').removeClass('text-primary').addClass('bg-primary')
        $('#pagination > *:nth-child(' + (step) + ') button span:first-child').addClass('bg-light').addClass('text-primary')
    }

    render() {
        let loading_output = ""
        if (this.state.loading)
            loading_output = <div className="alert alert-success">
                <div className="spinner-grow spinner-grow-sm" role="status">
                    <span className="visually-hidden">...</span>
                </div>
                <span className="ms-2">Loading parts</span>
            </div>
        let next_disabled = 'disabled'
        let next_text = "End assembly to continue"
        if (this.state.complete) {
            if (this.state.name === "") {
                next_text = "Set a name to continue"
            } else {
                next_disabled = ''
                next_text = "Create plasmid"
            }
        }

        let get_params = "n=" + this.state.name
        if (this.state.description)
            get_params += "&d=" + this.state.description
        if (this.state.receiver)
            get_params += "&b=" + this.state.receiver.i
        if (this.state.inserts.length) {
            let insert_names = []
            this.state.inserts.forEach((insert) => {
                insert_names.push(insert.i)
            })
            get_params += "&i=" + insert_names.join("+")
        }
        const enzymes_output = []
        AssemblyStandards[this.state.assembly_standard].enzymes.forEach((enzyme_data) => {
            let active = false
            if (enzyme_data.value === this.state.enzyme) active = true
            enzymes_output.push(<EnzymeRender name={enzyme_data.name} value={enzyme_data.value} active={active}
                                              setEnzymeHandler={this.setEnzymeHandler}/>)
        })
        const assembly_standards_options = []
        for (const [key, value] of Object.entries(AssemblyStandards)) {
            let selected = false
            if (this.state.assembly_standard === key) selected = true
            assembly_standards_options.push(<option selected={selected} value={key}>{value.name}</option>)
        }
        let api_error_output = ""
        if (this.state.apiError)
            api_error_output = <div className="alert alert-danger">{this.state.apiError}</div>

        // Part selector

        let receiver_output = <div>
            <h3>Receiver</h3>
            <p className="alert alert-info">Set name and select an enzyme to show options</p>
        </div>
        let inserts_output = <div>
            <h3>Inserts</h3>
            <p className="alert alert-info">Set a receiver to show options</p>
        </div>
        if (this.state.enzyme && this.state.name && this.state.parts && !this.state.loading) {
            receiver_output = []
            let receivers = []
            receiver_output.push(<h3>Receiver</h3>)
            receiver_output.push(<TableFilter text={this.state.receiverFilter}
                                              filterClearHandler={this.filterReceiverClearHandler}
                                              filterHandler={this.filterReceiverHandler}/>)
            this.state.parts.forEach((part) => {
                if (part.t === 1 && part.n.toLowerCase().includes(this.state.receiverFilter.toLowerCase())) {
                    let active = false
                    if (part === this.state.receiver) active = true
                    receivers.push(<PartRender assembly_standard={this.state.assembly_standard} part={part}
                                               active={active} config={this.state} button_text="Set"
                                               partHandler={this.setReceiverHandler}/>)
                }
            })
            if (receivers.length)
                receiver_output.push(<div className="aw_part_table mb-4">
                    <div className="alert alert-light">
                        <table className="table small">
                            <thead>
                            <tr>
                                <th scope="col">Name</th>
                                <th scope="col">Level - Type</th>
                                <th scope="col">Description</th>
                                <th scope="col">Set</th>
                            </tr>
                            </thead>
                            <tbody>{receivers}</tbody>
                        </table>
                    </div>
                </div>)
            else
                receiver_output.push(<div className="alert alert-danger">No receivers found</div>)

            inserts_output = []
            let inserts = []
            inserts_output.push(<h3>Inserts</h3>)

            if (!this.state.complete) {
                inserts_output.push(<TableFilter text={this.state.insertFilter}
                                                 filterClearHandler={this.filteriInsertClearHandler}
                                                 filterHandler={this.filteriInsertHandler}/>)
                let last_part_added = null
                if(this.state.receiver)
                    last_part_added = this.state.receiver
                if (this.state.inserts.length) {
                    last_part_added = this.state.inserts[this.state.inserts.length - 1]
                }
                this.state.parts.forEach((part) => {
                    if (part.t === 0 && part.n.toLowerCase().includes(this.state.insertFilter.toLowerCase())) {
                        if(last_part_added){
                            if(last_part_added.o3 !== part.o5){
                                return;
                            }
                        }
                        inserts.push(<PartRender assembly_standard={this.state.assembly_standard} part={part}
                                                 button_text="Add"
                                                 partHandler={this.addInsertHandler}/>)
                    }
                })
                if (inserts.length) {
                    inserts_output.push(<div className="aw_part_table mb-4">
                        <div className="alert alert-light">
                            <table className="table">
                                <thead>
                                <tr>
                                    <th scope="col">Name</th>
                                    <th scope="col">Level - Type</th>
                                    <th scope="col">Description</th>
                                    <th scope="col">Add</th>
                                </tr>
                                </thead>
                                <tbody>{inserts}</tbody>
                            </table>
                            <p className="alert alert-info fs-6">Only compatible inserts are listed</p>
                        </div>
                    </div>)
                } else
                    inserts_output.push(<div className="alert alert-danger">No inserts found</div>)

            } else {
                inserts_output.push(<div className="alert alert-success">Assembly complete</div>)
            }
        }

        let receiver_insert_output = []
        receiver_insert_output.push(<div className="row collapse">
            <Pagination step={2} position="t" stepHandler={this.stepHandler} create_url={window.create_url}></Pagination>
            {receiver_output}
            <Pagination step={2} position="b" stepHandler={this.stepHandler} create_url={window.create_url}></Pagination>
        </div>)
        receiver_insert_output.push(<div className="row collapse">
            <Pagination step={3} position="t" stepHandler={this.stepHandler} create_url={window.create_url} next_text={next_text}
                        next_disabled={next_disabled}></Pagination>
            {inserts_output}
            <Pagination step={3} position="b" stepHandler={this.stepHandler} create_url={window.create_url} next_text={next_text}
                        next_disabled={next_disabled}></Pagination>
        </div>)

        return <div id="assembly_wizard">
            <div className="container">
                <div className="row">
                    <div id="assembly_wizard-main">
                        {api_error_output}
                        <div className="row">
                            <div className="col-8">
                                <h2>Setup</h2>
                                <nav aria-label="Page navigation">
                                    <ul id="pagination" className="pagination">
                                        <li className="page-item active" onClick={() => {
                                            this.stepHandler(1)
                                        }}><button className="page-link"><span
                                            className="badge bg-light text-primary">Step 1</span><span className="ms-1">Name & enzyme</span></button>
                                        </li>
                                        <li className="page-item" onClick={() => {
                                            this.stepHandler(2)
                                        }}><button className="page-link"><span
                                            className="badge bg-primary">Step 2</span><span
                                            className="ms-1">Receiver</span></button>
                                        </li>
                                        <li className="page-item" onClick={() => {
                                            this.stepHandler(3)
                                        }}><button className="page-link"><span
                                            className="badge bg-primary">Step 3</span><span
                                            className="ms-1">Inserts</span></button></li>
                                    </ul>
                                </nav>
                                <div id="steps" className="border alert alert-light mb-4">
                                    <div className="collapse show">
                                        <div className="row">
                                            <Pagination step={1} position="t" stepHandler={this.stepHandler}
                                                        create_url={window.create_url}></Pagination>
                                            <div className="col-4">
                                                <h4>Name</h4>
                                                <p><input className="form-control" type="text"
                                                          placeholder="Plasmid name"
                                                          onChange={this.setName}/></p>
                                                <h4>Description</h4>
                                                <p><input className="form-control" type="text"
                                                          placeholder="Plasmid description"
                                                          onChange={this.setDescription}/></p>
                                            </div>
                                            <div className="col-4">
                                                <h4>Assembly Standard</h4>
                                                <div>
                                                    <select onChange={this.setStandardHandler} className="form-select">
                                                        {assembly_standards_options}
                                                    </select>
                                                </div>
                                            </div>
                                            <div className="col-4">
                                                <h4>Level / Enzyme</h4>
                                                <div>
                                                    <p>{enzymes_output}</p>
                                                </div>
                                            </div>
                                            <Pagination step={1} position="b" stepHandler={this.stepHandler}
                                                        create_url={window.create_url}></Pagination>
                                        </div>
                                    </div>
                                    {receiver_insert_output}
                                </div>
                                {loading_output}
                                <h2>Map</h2>
                                <PlasmidMap config={this.state} name={this.state.name}/>
                            </div>
                            <div className="col-4">
                                <h2>Assembly</h2>
                                <AssemblyResult create_url={window.create_url} next_text={next_text}
                                                next_disabled={next_disabled} config={this.state}
                                                removeLastPartHandler={this.removeLastPartHandler}/>
                                <form name={"form-wizard"} id={"form-wizard"} method="POST" action={window.next_url} className={"collapse"}>
                                    <input name="csrfmiddlewaretoken" value={this.state.csrf_token} type="hidden" />
                                    <input type={"hidden"} value={get_params} name={"params"} />
                                </form>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    }
}

export default App;
