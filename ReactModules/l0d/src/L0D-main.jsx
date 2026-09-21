import React, {Component} from 'react'
import Select from 'react-select'
import axios from 'axios';

import L0DPartsStandards from "./components/L0D-parts-standards"

import {
    clearSequence,
    getReverseComplementSequenceString,
    tmSeq
} from './components/Helpers'
import {
    aliasedEnzymesByName,
    getCutsitesFromSequence,
    getAminoAcidFromSequenceTriplet,
    aminoAcidToDegenerateDnaMap
} from 've-sequence-utils'

import Container from 'react-bootstrap/Container'
import Form from 'react-bootstrap/Form'
import FormControl from 'react-bootstrap/FormControl'
import FloatingLabel from 'react-bootstrap/FloatingLabel'
import OverlayTrigger from 'react-bootstrap/OverlayTrigger';
import Tooltip from 'react-bootstrap/Tooltip';


import $ from 'jquery';


class Highlight extends Component {
    render() {
        if (this.props.text) {
            let overlay_text = ""
            switch (this.props.type) {
                case 'eswen':
                    overlay_text = "Restriction site"
                    break
                case 'doh':
                    overlay_text = "Domestication overhang"
                    break
                case 'oh':
                    overlay_text = "Overhang"
                    break
                case 'seq':
                    overlay_text = "Sequence"
                    break
                case 'res':
                    overlay_text = "Internal restriction site"
                    break
                case 'res_dom':
                    overlay_text = "Internal restriction site (domesticated)"
                    break
                case 'prev_bases':
                    overlay_text = "Overhang complement bases for frame conservation"
                    break
                case 'stop':
                    overlay_text = "Stop codon"
                    break
                default:
                    break
            }
            if (this.props.extra)
                overlay_text += " " + this.props.extra
            return <OverlayTrigger key={"ot" + this.props.key} placement="top" overlay={
                <Tooltip id={"id" + this.props.key}>{overlay_text}</Tooltip>
            }>
                <span className={"hl hl-" + this.props.type}>{this.props.text}</span>
            </OverlayTrigger>
        }
        return <i></i>
    }
}

class Pagination extends React.Component {
    render() {
        let next_button = ""
        if (this.props.step < 5)
            next_button = <button className="btn btn-success float-end" onClick={() => {
                this.props.stepHandler(this.props.step + 1)
        }}>Next</button>

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
                {next_button}
                {prev_button}
            </div>
            {post_hr}
        </div>
    }
}


class InputOHS extends Component {
    render() {
        const standard = L0DPartsStandards[this.props.standard]
        return <Form.Select
            onChange={this.props.handler}
            value={this.props.cv}
            aria-label={standard.name + " OHs / " + this.props.oh + "'"}
            className="oh-select" ref={"ref" + this.props.oh}>
            {Object.keys(standard.ohs).map((key) => {
                return (
                    <option key={key}
                            value={key}>{standard.ohs[key].name + " - " + standard.ohs[key].oh + (standard.ohs[key].prev_bases ? ' [ ' + standard.ohs[key].prev_bases + ' ]' : '') + (standard.ohs[key].stop ? ' [ STOP ]' : '')}</option>
                )
            })}
            <option key="custom" value="custom">Custom</option>
        </Form.Select>
    }
}

class OHInput extends Component {
    render() {
        return <FloatingLabel controlId={"ohs" + this.props.oh}
                              label={this.props.standard + " OHs / " + this.props.oh + "'"}>
            <InputOHS
                standard={this.props.standard}
                oh={this.props.oh}
                handler={this.props.handler}
                cv={this.props.cv}
            />
        </FloatingLabel>
    }
}

class Primer extends Component {
    render() {
        return <tr>
            <th scope="row">{this.props.type} <i className="bi bi-clipboard copy_clipboard" data-cc={this.props.name + "-" + this.props.type}></i></th>
            <td>{this.props.seq.length}</td>
            <td><i className="bi bi-clipboard copy_clipboard" data-cc={this.props.seq}></i> {this.props.seq}</td>
        </tr>
    }
}

class PrimerDesign extends Component {
    render() {
        const fragment_name = this.props.name + "-F" + (this.props.idx + 1)
        let primers = []
        const finalSeq = this.props.eswen + this.props.extra5 + this.props.template + this.props.extra3 + getReverseComplementSequenceString(this.props.eswen)
        let method = "pcr"
        let method_text = "PCR"
        let tm_warning = ""
        if (finalSeq.length <= this.props.pcrMinLength) {
            method = "oann"
            method_text = "GBlock or Oligo Annealing"
        }
        if (method === "oann") {
            primers.push(<Primer seq={finalSeq} name={fragment_name} type={"F"}/>)
            primers.push(<Primer seq={getReverseComplementSequenceString(finalSeq)} name={fragment_name} type={"R"}/>)
        } else {
            let tmSeq_dir = tmSeq(this.props.template, this.props.pcrTm)
            let tmSeq_rc = tmSeq(getReverseComplementSequenceString(this.props.template), this.props.pcrTm)
            if (tmSeq_dir === "minTempNotReached" || tmSeq_dir === "minLengthError") {
                tmSeq_dir = this.props.template
                tmSeq_rc = getReverseComplementSequenceString(this.props.template)
                if (tmSeq_dir === "minTempNotReached")
                    tm_warning = <div className="alert alert-danger">Tm below min requirement</div>
                if (tmSeq_dir === "minLengthError")
                    tm_warning = <div className="alert alert-danger">Template length below min requirement (7 bp)</div>
            }
            const finalSeq_fwd = this.props.eswen + this.props.extra5 + tmSeq_dir
            primers.push(<Primer seq={finalSeq_fwd} name={fragment_name} type={"F"}/>)
            const finalSeq_rev = this.props.eswen + getReverseComplementSequenceString(this.props.extra3) + tmSeq_rc
            primers.push(<Primer seq={finalSeq_rev} name={fragment_name} type={"R"}/>)
            // primers.push(<Primer seq={getReverseComplementSequenceString(tmSeq_rc) + this.props.extra3 + getReverseComplementSequenceString(this.props.eswen) } name={fragment_name} type={"R (RC)"}/>)
        }
        return (<div className="alert alert-light border text-break">
            <h5>
                <span className="badge bg-secondary text-light">{fragment_name}</span>
                <span className="small fw-light ms-2">{method_text} / {finalSeq.length} bp</span>
            </h5>
            {tm_warning}
            <table className="table small">
                <thead>
                <tr>
                    <th scope="col">Type</th>
                    <th scope="col">Length</th>
                    <th scope="col">Seq</th>
                </tr>
                </thead>
                <tbody>
                    {primers}
                </tbody>
            </table>
        </div>)
    }
}

class ResEditor extends Component {
    constructor(props) {
        super(props)
        if (props.new_res === undefined)
            this.resInputHandle(this.props.res.idx, this.props.res.seq)
    }

    resInputHandle = (i, v) => {
        this.props.resInputHandle(i, v)
    }

    getNucleotidesFromIUPAC(letter) {
        switch (letter.toLowerCase()) {
            case 'a':
                return ['a']
            case 't':
                return ['t']
            case 'c':
                return ['c']
            case 'g':
                return ['g']
            case 'r':
                return ['g', 'a']
            case 'y':
                return ['t', 'c']
            case 'm':
                return ['a', 'c']
            case 'k':
                return ['g', 't']
            case 's':
                return ['g', 'c']
            case 'w':
                return ['a', 't']
            case 'h':
                return ['a', 'c', 't']
            case 'b':
                return ['g', 't', 'c']
            case 'v':
                return ['g', 'c', 'a']
            case 'd':
                return ['g', 'a', 't']
            case 'n':
                return ['g', 'a', 'c', 't']
            default:
                return letter
        }
    }

    generateVariations(triplet) {
        if (triplet.length === 1)
            return this.getNucleotidesFromIUPAC(triplet)

        const alternatives = []
        const firstLetterVariations = this.getNucleotidesFromIUPAC(triplet.charAt(0))
        const sub_alternatives = this.generateVariations(triplet.substring(1, triplet.length))

        firstLetterVariations.forEach((firstLetterVariation) => {
            sub_alternatives.forEach((sub_alternative) => {
                alternatives.push(firstLetterVariation + sub_alternative)
            })
        })
        return alternatives
    }

    getVariation(domesticate) {
        if (domesticate.length !== 3) {
            alert("Domestication site length != 3")
            return false
        }
        const triplet = aminoAcidToDegenerateDnaMap[getAminoAcidFromSequenceTriplet(domesticate).value.toLowerCase()]
        const results = this.generateVariations(triplet)
        let result = false
        results.forEach((v) => {
            if (v !== domesticate) {
                result = v
            }
        })
        return result
    }

    recommended_codon() {
        const new_res = this.props.new_res
        const res = this.props.res
        if (new_res && new_res.length >= 3) {
            let result = res.seq.substring(0, res.start % 3)
            const domesticate = res.seq.substring(res.start % 3, res.start % 3 + 3)
            const variation = this.getVariation(domesticate)
            if (variation)
                result += variation
            else {
                alert("No variation found")
                return res
            }
            result += res.seq.substring(res.start % 3 + 3, new_res.length)
            return result
        }
        alert("No Restriction enzyme site set")
        return ""
    }

    render() {
        if (this.props.res) {
            let equal_warwing = ""
            let res = this.props.res
            let new_res = this.props.new_res ? this.props.new_res : ""

            if (res.seq.toLowerCase() === new_res.toLowerCase())
                equal_warwing = <div className="alert alert-danger">Restriction site unchanged</div>
            let diferent_length_warwing = ""

            if (res.seq.length !== new_res.length && new_res.length > 0)
                diferent_length_warwing =
                    <div className="alert alert-warning">New restriction site has different length than original</div>

            let frame_string = ""
            if (new_res) {
                let start = new_res.substring(0, res.start % 3) + "-"
                if (start.length === 1)
                    start = ""
                frame_string = " (" + start + new_res.substring(res.start % 3, res.start % 3 + 3) + ")"
            }
            const recommended_codon = this.recommended_codon()

            return <div className="alert alert-info border">
                <h5>Restriction site # {res.idx + 1}</h5>
                <div className="row">
                    <div className="col-7">
                        <FloatingLabel controlId={"resInput" + res.idx} label="New restriction enzyme site">
                            <FormControl onChange={(e) => this.resInputHandle(res.idx, e.target.value)}
                                         value={new_res} aria-label="New restriction enzyme site"/>
                        </FloatingLabel>
                        <div className="mb-2">
                            <button type="button" className="btn btn-secondary me-2"
                                    onClick={(e) => this.resInputHandle(res.idx, res.seq)}>Reset
                            </button>
                            <button type="button" className="btn btn-success me-2"
                                    onClick={(e) => this.resInputHandle(res.idx, recommended_codon)}>Recommended
                            </button>
                        </div>
                    </div>
                    <div className="col-5">
                        <div className="mb-2 small">
                            <span><strong>Enzyme</strong>: {res.enzyme}</span><br/>
                            <span><strong>Position</strong>: {(res.start + 1) + " to " + (res.end + 1)}</span><br/>
                            <span><strong>Frame position</strong>: {((res.start % 3) + 1) + frame_string}</span><br/>
                        </div>
                    </div>
                </div>
                {equal_warwing}
                {diferent_length_warwing}
            </div>
        } else {
            return <div></div>
        }
    }
}

class Fragments extends Component {
    constructor(props) {
        super(props)
        let new_res = []
        props.res.forEach((re) => {
            new_res.push(re.seq)
        })
        this.state = {
            new_res: new_res,
            apiCall_result: []
        }
    }


    componentWillReceiveProps(nextProps) {
        let new_res = []
        nextProps.res.forEach((re) => {
            new_res.push(re.seq)
        })
        this.setState({
            new_res: new_res,
            apiCall_result: []
        })
    }

    commonStart(seq1, seq2) {
        let finalSeq = ""
        for (var i = 0; i < Math.min(seq1.length, seq2.length); i++) {
            if (seq1.charAt(i).toLowerCase() === seq2.charAt(i).toLowerCase()) {
                finalSeq += seq1.charAt(i)
            } else
                break
        }
        return finalSeq
    }

    firstCommonIndexAfterDiff(ref, query) {
        let firstIndex = null
        let diffFound = false
        if (ref.length === query.length)
            for (var i = 0; i < query.length; i++) {
                if (ref.charAt(i).toLowerCase() === query.charAt(i).toLowerCase()){
                    if (diffFound && ref.substring(i, query.length) === query.substring(i, query.length)) {
                        firstIndex = i
                        break
                    }
                } else
                        diffFound = true
            }
        return firstIndex
    }

    lastCommonIndex(ref, query) {
        let lastIndex = 0
        for (var i = 0; i < Math.min(ref.length, query.length); i++) {
            if (ref.charAt(i).toLowerCase() === query.charAt(i).toLowerCase())
                lastIndex = i
            else
                break
        }
        return lastIndex
    }

    resParts(ref, query) {
        const lastCommonIndex = this.lastCommonIndex(ref, query)
        const firstCommonIndexAfterDiff = this.firstCommonIndexAfterDiff(ref, query)

        // console.log("lastCommonIndex: " + lastCommonIndex)
        // console.log("firstCommonIndexAfterDiff: " + firstCommonIndexAfterDiff)

        let start = ref.substring(0, lastCommonIndex + 1)
        let diff = ""
        let end = ""
        if (firstCommonIndexAfterDiff) {
            diff = ref.substring(lastCommonIndex + 1, firstCommonIndexAfterDiff)
            end = ref.substring(firstCommonIndexAfterDiff, ref.length)
        } else {
            diff = ref.substring(lastCommonIndex + 1, ref.length)            
        }

        return {
            'start': start,
            'diff': diff,
            'end': end,
        }
    }

    resInputHandle = (i, v) => {
        let new_res = this.state.new_res
        new_res[i] = clearSequence(v)
        this.setState({
            new_res: new_res,
            apiCall_result: []
        })
    }

    formatOHEnd(seq, len) {
        return seq.substring(0, seq.length - len).toLowerCase() + seq.substring(seq.length - len, seq.length).toUpperCase()
    }

    formatOH(seq, len) {
        return seq.substring(0, len).toUpperCase() + seq.substring(len, seq.length).toLowerCase()
    }

    stepHandler = (step) => {
        this.props.stepHandler(step)
    }

    apiCall(ohs) {
        this.setState({
            apiCall_result: {
                'wait': true
              }
        })    
        // const url = 'http://192.168.9.61:8000/inventory/api/fidelity_calc/sapi/' + ohs.join('-')
        const url = '/inventory/api/fidelity_calc/sapi/' + ohs.join('-')
        axios.get(url)
            .then((response) => {
                var replaced = response.data.ligation_frequency_matrix_html.replaceAll('container', 'lfm_container')
                replaced = replaced.replaceAll('row', 'lfm_row')
                replaced = replaced.replaceAll('cell', 'lfm_cell')
                this.setState({
                    apiCall_result: {
                        'fidelity': response.data.fidelity,
                        'render': replaced
                      }
                })
            })
            .catch((error) => {
                this.setState({
                    apiCall_result: {
                        'error': error
                      }
                })              
            })
    }

    render() {
        let amplicon_output = []
        let fragments_output = []
        let final_sequence = []
        let final_sequence_txt = ""
        let domestication_output = []
        domestication_output.push(<h3>Restriction sites</h3>)
        let domesticaions_to_make = false
        let fragments = []
        let extra5 = ""
        let extra3 = ""
        let minOhLengthAlert = ""
        let ohs = []
        let domestication_pending = false
        const stop_codon = 'tga'
        const oh_length = Math.abs(this.props.the_re.topSnipOffset - this.props.the_re.bottomSnipOffset)

        ohs.push(this.props.l0_receiver.ohs.oh5)

        final_sequence.push(<Highlight text={this.props.eswen} type="eswen" key="1"
                                       extra={" (" + this.props.the_re.name + ")"}/>)
        final_sequence_txt += this.props.eswen
        final_sequence.push(<Highlight text={this.props.l0_receiver.ohs.oh5} type="doh" key="2"/>)
        final_sequence_txt += this.props.l0_receiver.ohs.oh5
        extra5 += this.props.l0_receiver.ohs.oh5.toUpperCase()

        final_sequence.push(<Highlight text={this.props.oh5.oh} extra={this.props.oh5.name} type="oh" key="3"/>)
        final_sequence_txt += this.props.oh5.oh
        extra5 += this.props.oh5.oh.toUpperCase()

        if (this.props.res.length) {
            let from = 0
            this.props.res.forEach((re, idx) => {
                final_sequence.push(<Highlight text={this.props.seq.substring(from, re.start)} type="seq"
                                               extra={" fragment " + (idx + 1)} key={"s-" + idx}/>)
                final_sequence_txt += this.props.seq.substring(from, re.start)
                let template = this.props.seq.substring(from, re.start)

                let type = "res"
                let seq = re.seq
                if (this.state.new_res[idx] && this.state.new_res[idx] !== re.seq) {
                    type = "res_dom"
                    seq = this.state.new_res[idx]
                }
                final_sequence.push(<Highlight text={seq} type={type} extra={" # " + (idx + 1)} key={"r-" + idx}/>)
                final_sequence_txt += seq
                template += this.commonStart(this.state.new_res[idx], re.seq)

                if (idx) {
                    // not the first
                    const resParts = this.resParts(this.state.new_res[idx - 1], this.props.res[idx - 1].seq)
                    extra5 = this.formatOH(resParts.diff, oh_length)
                    template = this.formatOH(resParts.end + template, oh_length - extra5.length)
                    ohs.push((resParts.diff + template).substring(0, oh_length).toUpperCase())
                }

                const resParts = this.resParts(this.state.new_res[idx], re.seq)
                extra3 = this.formatOH((resParts.diff + resParts.end + this.props.seq.substring(re.end + 1, this.props.seq.length)).substring(0, oh_length), oh_length)
                if (idx !== this.props.res.length - 1)
                    template = this.formatOHEnd(template, oh_length - extra3.length)
                if (extra3.length < oh_length)
                    minOhLengthAlert = <div className="alert alert-danger">OH length under optimal length</div>

                fragments.push({
                    idx: idx,
                    extra5: extra5,
                    template: template,
                    extra3: extra3,
                })
                domestication_output.push(<ResEditor resInputHandle={this.resInputHandle}
                                                     new_res={this.state.new_res[idx]} res={re}/>)
                if (re.seq.toLowerCase() === this.state.new_res[idx].toLowerCase())
                    domestication_pending = true
                domestication_output.push(minOhLengthAlert)
                domesticaions_to_make = true
                from = re.end + 1

                if (idx === this.props.res.length - 1) {
                    // the last
                    final_sequence.push(<Highlight text={this.props.seq.substring(re.end + 1, this.props.seq.length)}
                                                   type="seq"
                                                   extra={" fragment " + (idx + 2)} key={"s-" + (idx + 1)}/>)
                    final_sequence_txt += this.props.seq.substring(re.end + 1, this.props.seq.length)
                    const resParts = this.resParts(this.state.new_res[idx], re.seq)
                    extra5 = this.formatOH(resParts.diff, oh_length)
                    template = this.formatOH(resParts.end + this.props.seq.substring(from, this.props.seq.length), oh_length - extra5.length)
                    ohs.push((resParts.diff + template).substring(0, oh_length).toUpperCase())
                    fragments.push({
                        idx: idx + 1,
                        extra5: extra5,
                        template: template,
                        extra3: "",
                    })
                }
            })
        } else {
            final_sequence.push(<Highlight text={this.props.seq} type="seq" key="s-1"/>)
            final_sequence_txt += this.props.seq
            fragments.push({
                idx: 0,
                extra5: extra5,
                template: this.props.seq.toLowerCase(),
                extra3: "",
            })
        }

        final_sequence.push(<Highlight text={this.props.oh3.prev_bases ? this.props.oh3.prev_bases : ''} type="prev_bases" key="6"/>)
        final_sequence_txt += this.props.oh3.prev_bases ? this.props.oh3.prev_bases : ''
        extra3 = this.props.oh3.prev_bases ? this.props.oh3.prev_bases : ''
        final_sequence.push(<Highlight text={this.props.oh3.stop ? stop_codon : ''} type="stop" key="7"/>)
        final_sequence_txt += this.props.oh3.stop ? stop_codon : ''
        extra3 += this.props.oh3.stop ? stop_codon : ''
        final_sequence.push(<Highlight text={this.props.oh3.oh} extra={this.props.oh3.name} type="oh" key="8"/>)
        final_sequence_txt += this.props.oh3.oh
        extra3 += this.props.oh3.oh.toUpperCase()
        final_sequence.push(<Highlight text={this.props.l0_receiver.ohs.oh3} type="doh" key="9"/>)
        final_sequence_txt += this.props.l0_receiver.ohs.oh3
        extra3 += this.props.l0_receiver.ohs.oh3.toUpperCase()
        
        ohs.push(this.props.l0_receiver.ohs.oh3)
        
        final_sequence.push(<Highlight text={getReverseComplementSequenceString(this.props.eswen)} type="eswen"
                                       key="10" extra={" (RevComp) (" + this.props.the_re.name + ")"}/>)
        final_sequence_txt += getReverseComplementSequenceString(this.props.eswen)
        fragments[fragments.length - 1].extra3 = extra3

        fragments.forEach((fragment) => {
            fragments_output.push(<PrimerDesign name={this.props.name} pcrMinLength={this.props.pcrMinLength}
                                                pcrTm={this.props.pcrTm} idx={fragment.idx} extra5={fragment.extra5}
                                                extra3={fragment.extra3} template={fragment.template}
                                                eswen={this.props.eswen}/>)
        })

        if (!domesticaions_to_make)
            domestication_output.push(<div className="alert alert-info">No domestication sites</div>)

        domestication_output.push(<h3>Ligation overhangs<i className="bi bi-clipboard copy_clipboard ms-2" data-cc={ohs.join(",")}></i></h3>)
        
        let ligation_fidelity_output = []
        ligation_fidelity_output.push(<div className="alert alert-warning">Overhangs: {ohs.join(" + ")}</div>)
        if(this.state.apiCall_result['wait']){
            ligation_fidelity_output.push(<button className="btn btn-success" onClick={() => {this.apiCall()}}>Wait <span class="spinner-grow spinner-grow-sm text-light" role="status"><span class="visually-hidden">Loading...</span></span></button>)
        } else {
            if(this.state.apiCall_result['fidelity']){
                ligation_fidelity_output.push(<div className="alert alert-success">Ligation Fidelity: {this.state.apiCall_result['fidelity']}%</div>)
                ligation_fidelity_output.push(<div className="alert alert-light">
                        <h4>Ligation fidelity matrix</h4>
                        <div id="fidelity_result" dangerouslySetInnerHTML={{ __html: this.state.apiCall_result['render'] }} />
                    </div>)  
                ligation_fidelity_output.push(<div className="alert alert-info">Conditions: SapI / 37-16 cycling. Using <a href="https://ligasefidelity.neb.com/viewset/run.cgi" target="_blank" rel="noreferrer">NEBridge Ligase Fidelity Viewer <i className="bi bi-box-arrow-up-right"></i></a></div>)
            } else {
                if(this.state.apiCall_result['error']){
                    ligation_fidelity_output.push(<div className="alert alert-danger">Error. Try again. [{this.state.apiCall_result['error']}]</div>)
                }
                ligation_fidelity_output.push(<button className="btn btn-success" onClick={() => {this.apiCall(ohs)}}>Estimate fidelity</button>)
            }
        }
        domestication_output.push(<div className="alert alert-light border text-break">
            <div className="mb-2">{ligation_fidelity_output}</div>
        </div>)
        domestication_output.push(<Pagination step={3} position="b" stepHandler={this.stepHandler}></Pagination>)



        let final_sequence_no_eswen_ohs_txt = final_sequence_txt.substring(this.props.eswen.length + this.props.l0_receiver.ohs.oh5.length + this.props.oh5.oh.length, final_sequence_txt.length - this.props.eswen.length - this.props.l0_receiver.ohs.oh3.length - this.props.oh3.oh.length)
        if (domestication_pending){
            amplicon_output.unshift(<div className="alert alert-danger">Restriction site domestication(s) pending</div>)
            fragments_output.unshift(<div className="alert alert-danger">Restriction site domestication(s) pending</div>)
        }
        amplicon_output.unshift(<h3>{this.props.name} <i className="bi bi-clipboard copy_clipboard ms-2" data-cc={final_sequence_txt}></i> <span className="badge text-dark">
            <form method="POST" action="/inventory/plasmid/create/l0d" target="_blank" class="default-style">
                <input type="hidden" name="oh5-name" value={this.props.oh5.name} />
                <input type="hidden" name="oh5-oh" value={this.props.oh5.oh} />
                <input type="hidden" name="oh3-name" value={this.props.oh3.name} />
                <input type="hidden" name="oh3-prev_bases" value={this.props.oh3.prev_bases ? this.props.oh3.prev_bases : ''} />
                <input type="hidden" name="oh3-stop" value={this.props.oh3.stop ? stop_codon : ''} />
                <input type="hidden" name="oh3-oh" value={this.props.oh3.oh} />
                <input type="hidden" name="name" value={this.props.name} />
                <input name="csrfmiddlewaretoken" value={window.csrf_token} type="hidden" />
                <input type="hidden" name="seq" value={final_sequence_no_eswen_ohs_txt} />
                <button className="btn btn-primary btn-sm" type="submit">Create plasmid <i className="bi bi-box-arrow-up-right"></i></button>
            </form>
        </span></h3>)

        amplicon_output.push(<div className="alert alert-light border text-break">{final_sequence}</div>)
        amplicon_output.push(<Pagination step={4} position="b" stepHandler={this.stepHandler}></Pagination>)

        fragments_output.push(<Pagination step={5} position="b" stepHandler={this.stepHandler}></Pagination>)

        return [
            <div className="collapse">{domestication_output}</div>,
            <div className="collapse text-break">{amplicon_output}</div>,
            <div className="collapse">{fragments_output}</div>,
            ]
    }
}

const default_standard = 'ytk'
const initialState = {
    sequence: '',
    standard: default_standard,
    custom_oh5: '',
    custom_oh3: '',
    oh5: L0DPartsStandards[default_standard].default[5],
    oh3: L0DPartsStandards[default_standard].default[3],
    pcrTm: 60,
    pcrMinLength: 80,
    name: 'Part',
    base_pairs_from_end: 'aa',
    enzymes: L0DPartsStandards[default_standard].domestication_enzymes,
}

class L0D extends Component {
    constructor(props) {
        super(props)
        this.state = initialState
    }

    resetState = () => {
        this.setState(initialState);
    }

    basePairsFromEndInputChangeHandle = (event) => {
        this.setState({
            base_pairs_from_end: event.target.value,
        })
    }


    domEnzymesInputChangeHandle = (event) => {
        let newEnzymes = []
        event.forEach((enzyme) => {
            newEnzymes.push(enzyme.value)
        })
        this.setState({
            enzymes: newEnzymes,
        })
    }

    pcrMinLengthInputChangeHandle = (event) => {
        this.setState({
            pcrMinLength: event.target.value,
        })
    }

    pcrTmInputChangeHandle = (event) => {
        this.setState({
            pcrTm: event.target.value,
        })
    }

    nameInputChangeHandle = (event) => {
        this.setState({
            name: event.target.value,
        })
    }

    sequenceInputChangeHandle = (event) => {
        this.setState({
            sequence: clearSequence(event.target.value),
        })
    }

    partStandardChangeHandle = (event) => {
        const standard = L0DPartsStandards[event.target.value]
        if (standard)
            this.setState({
                standard: event.target.value,
                oh5: standard.default[5],
                oh3: standard.default[3],
                enzymes: standard.domestication_enzymes
            })
    }

    OH5InputChangeHandle = (event) => {
        this.setState({oh5: event.target.value})
    }

    OH3InputChangeHandle = (event) => {
        this.setState({oh3: event.target.value})
    }

    customOh5InputHandle = (event) => {
        this.setState({custom_oh5: clearSequence(event.target.value)})
    }

    customOh3InputHandle = (event) => {
        this.setState({custom_oh3: clearSequence(event.target.value)})
    }

    partStandardInputItems = Object.keys(L0DPartsStandards).map(function (key) {
        return <option key={key} value={key}>{L0DPartsStandards[key].name}</option>
    })

    getEnzymesSelect = (enzymes) => {
        let domesticationEnzymesOutput = []
        enzymes.forEach((enzyme) => {
            domesticationEnzymesOutput.push({value: enzyme, label: aliasedEnzymesByName[enzyme].name})
        })
        return domesticationEnzymesOutput
    }

    findInFrame(targets, query) {
        for (let t = 0; t < targets.length; t++) {
            for (let q = 0; q < query.length; q = q + 3) {
                if (targets[t].toLowerCase() === query.substring(q, q + 3).toLowerCase()) {
                    return q
                }
            }
        }
        return false
    }

    stepHandler(step) {
        $('#steps').children().removeClass('show')
        $('#steps > *:nth-child(' + (step) + ')').addClass('show')

        let li = $('#pagination').children()
        let li_activo = $('#pagination > *:nth-child(' + (step) + ')')

        li.removeClass('active')
        li_activo.addClass('active')

        li.each(function(){
            $(this).find('span:first-child').removeClass('bg-light').removeClass('text-primary').addClass('bg-'+$(this).attr('data-color'))
            $(this).find('button span:last-child').addClass('text-'+$(this).attr('data-color'))
        })
        li_activo.find('button span:first-child').addClass('bg-light').addClass('text-primary')
        li_activo.find('button span:last-child').removeClass('text-'+li_activo.attr('data-color'))
    }

    componentDidUpdate() {
        window.onReady()
    }

    componentDidMount() {
        // this.stepHandler(3)
    }

    render() {
        const the_standard = L0DPartsStandards[this.state.standard]
        const domesticationEnzymesOptions = []
        Object.keys(aliasedEnzymesByName).forEach(function (k, v) {
            if (aliasedEnzymesByName[k].isType2S)
                domesticationEnzymesOptions.push({value: k, label: aliasedEnzymesByName[k].name})
        })

        let custom_oh5_input
        if (this.state.oh5 === "custom") {
            custom_oh5_input = <FloatingLabel controlId="customOh5Input" label="Custom OH 5'">
                <FormControl onChange={this.customOh5InputHandle} value={this.state.custom_oh5}
                             aria-label="Custom OH 5"/>
            </FloatingLabel>
        }

        let custom_oh3_input
        if (this.state.oh3 === "custom") {
            custom_oh3_input = <FloatingLabel controlId="customOh3Input" label="Custom OH 3'">
                <FormControl onChange={this.customOh3InputHandle} value={this.state.custom_oh3}
                             aria-label="Custom OH 3"/>
            </FloatingLabel>
        }

        const sequenceInput = this.state.sequence.toLowerCase()
        
        let output = []

        if (!sequenceInput) {
            output = [
                <div className="collapse">
                <div className="alert alert-info">Set a sequence to continue</div>
                <Pagination step={3} position="b" stepHandler={this.stepHandler}></Pagination>
                </div>,
                <div className="collapse">
                <div className="alert alert-info">Set a sequence to continue</div>
                <Pagination step={4} position="b" stepHandler={this.stepHandler}></Pagination>
                </div>,
                <div className="collapse">
                <div className="alert alert-info">Set a sequence to continue</div>
                <Pagination step={5} position="b" stepHandler={this.stepHandler}></Pagination>
                </div>,
                ]
        } else {
            let RES = []
            this.state.enzymes.forEach((enzyme_name) => {
                RES.push(aliasedEnzymesByName[enzyme_name])
            })
            const cutSites = getCutsitesFromSequence(sequenceInput, true, RES)

            let restrictionSites = []

            if (Object.keys(cutSites).length) {
                Object.entries(cutSites).forEach(([key, enzymeCuts]) => {
                    enzymeCuts.forEach((enzymeCut, idx) => {
                        const start = enzymeCut.recognitionSiteRange['start']
                        const end = enzymeCut.recognitionSiteRange['end']
                        if (start < end)
                            restrictionSites.push({
                                seq: sequenceInput.substring(enzymeCut.recognitionSiteRange['start'], enzymeCut.recognitionSiteRange['end'] + 1),
                                start: start,
                                end: end,
                                enzyme: enzymeCut.name
                            })
                    })
                    restrictionSites.sort(function (a, b) {
                        return a.start - b.start
                    })
                    restrictionSites.forEach((restrictionSite, idx) => {
                        restrictionSite.idx = idx
                    })
                })
            }

            const the_re = aliasedEnzymesByName[the_standard.enzyme]
            const enzymeSiteWithExtraNucl = this.state.base_pairs_from_end + the_re.site.toUpperCase() + the_standard.bases_upto_snip

            let oh5 = {
                name: 'Custom',
                oh: this.state.custom_oh5
            }
            if (this.state.oh5 !== "custom") {
                oh5 = the_standard.ohs[this.state.oh5]
            }
            let oh3 = {
                name: 'Custom',
                oh: this.state.custom_oh3
            }
            if (this.state.oh3 !== "custom") {
                oh3 = the_standard.ohs[this.state.oh3]
            }

            output = <Fragments name={this.state.name} pcrMinLength={this.state.pcrMinLength} pcrTm={this.state.pcrTm}
                                the_re={the_re} res={restrictionSites} oh5={oh5} oh3={oh3}
                                l0_receiver={the_standard.l0_receiver} eswen={enzymeSiteWithExtraNucl} seq={sequenceInput}
                                stepHandler={this.stepHandler}/>

        }

        let atg_output = <div className="alert alert-alert border">No ATG found in frame 1</div>
        let stop_output = <div className="alert alert-light border">No STOP codon found in frame 1</div>

        let atg_found = this.findInFrame(['atg'], sequenceInput)
        let stop_found = this.findInFrame(['taa', 'tga', 'tag'], sequenceInput)

        if (atg_found !== false)
            atg_output = <div className="alert alert-success">ATG found in frame 1, position {atg_found + 1}</div>

        if (stop_found !== false)
            stop_output =
                <div className="alert alert-danger">STOP codon found in frame 1, position {stop_found + 1}</div>

        return (
            <Container>

                <div className="flow_root">
                    <nav aria-label="Page navigation" className="float-start">
                        <ul id="pagination" className="pagination">
                            <li className="page-item active" data-color="primary" onClick={() => {
                                this.stepHandler(1)
                            }}><button className="page-link"><span
                                className="badge bg-light text-primary">Step 1</span><span
                                className="ms-1">Standard, enzymes & params</span></button>
                            </li>
                            <li className="page-item" data-color="primary" onClick={() => {
                                this.stepHandler(2)
                            }}><button className="page-link"><span
                                className="badge bg-primary">Step 2</span><span
                                className="ms-1">Sequence</span></button>
                            </li>
                            <li className="page-item" data-color="primary" onClick={() => {
                                this.stepHandler(3)
                            }}><button className="page-link"><span
                                className="badge bg-primary">Step 3</span><span
                                className="ms-1">Domestication</span></button></li>
                            <li className="page-item" data-color="success" onClick={() => {
                                this.stepHandler(4)
                            }}><button className="page-link"><span
                                className="badge bg-success">Result 1</span><span
                                className="ms-1 text-success fw-bold">Final sequence</span></button></li>
                            <li className="page-item" data-color="success" onClick={() => {
                                this.stepHandler(5)
                            }}><button className="page-link"><span
                                className="badge bg-success">Result 2</span><span
                                className="ms-1 text-success fw-bold">Oligos</span></button></li>
                        </ul>
                    </nav>
                    <button className="btn btn-secondary me-2 float-end" onClick={() => {
                            this.resetState()
                    }}><i class="bi bi-arrow-clockwise"></i></button>
                </div>
                <div id="steps" className="border alert alert-light mb-4">
                    <div className="collapse show">
                        <h3>Standard</h3>
                        <FloatingLabel controlId="partStandardInput" label="Assembly Standard">
                            <Form.Select onChange={this.partStandardChangeHandle} aria-label="Part standard input">
                                {this.partStandardInputItems}
                            </Form.Select>
                        </FloatingLabel>
                        <h3 className="mt-2">Position</h3>
                        <OHInput standard={this.state.standard} oh="5" cv={this.state.oh5}
                                 handler={this.OH5InputChangeHandle} />
                        {custom_oh5_input}
                        <OHInput standard={this.state.standard} oh="3" cv={this.state.oh3}
                                 handler={this.OH3InputChangeHandle} />
                        {custom_oh3_input}
                        <h3 className="mt-2">Enzymes</h3>
                        <Form.Text className="text-muted">
                            Domestication enzymes
                        </Form.Text>
                        <Select options={domesticationEnzymesOptions} value={this.getEnzymesSelect(this.state.enzymes)}
                                isMulti className="basic-multi-select mb-2" classNamePrefix="select"
                                onChange={this.domEnzymesInputChangeHandle}/>
                        <h3 className="mt-2">Parameters</h3>
                        <div>
                            <FloatingLabel controlId="pcrMinLengthInput"
                                           label="Length under which part is prepared by GBlock or oligo annealing instead of PCR">
                                <FormControl onChange={this.pcrMinLengthInputChangeHandle}
                                             value={this.state.pcrMinLength} aria-label="pcr/oann min length"/>
                            </FloatingLabel>
                            <FloatingLabel controlId="pcrTmInput" label="Tm PCR primers °C">
                                <FormControl onChange={this.pcrTmInputChangeHandle}
                                             value={this.state.pcrTm} aria-label="PCR Tm primers"/>
                            </FloatingLabel>
                            <FloatingLabel controlId="basePairsFromEnd" label="Base pairs from end">
                                <FormControl onChange={this.basePairsFromEndInputChangeHandle}
                                             value={this.state.base_pairs_from_end} aria-label="Base pairs from end"/>
                            </FloatingLabel>
                        </div>
                        <Pagination step={1} position="b" stepHandler={this.stepHandler}></Pagination>
                    </div>
                    <div className="collapse">
                        <h3>Part name</h3>
                        <FloatingLabel controlId="nameInput" label="Name Input">
                            <FormControl onChange={this.nameInputChangeHandle}
                                         value={this.state.name} aria-label="Name Input"/>
                        </FloatingLabel>
                        <h3 className="mt-2">Sequence ({this.state.sequence.length} bp)</h3>
                        <FloatingLabel controlId="sequenceInput" label="Sequence Input">
                            <FormControl onChange={this.sequenceInputChangeHandle} as="textarea"
                                         value={this.state.sequence} aria-label="Sequence Input"/>
                            <Form.Text className="text-muted">
                                Non ATGC characters are automatically removed.
                            </Form.Text>
                        </FloatingLabel>
                        {atg_output}
                        {stop_output}
                        <Pagination step={2} position="b" stepHandler={this.stepHandler}></Pagination>
                    </div>
                    {output}
                </div>
            </Container>
        )
    }
}

export default L0D
