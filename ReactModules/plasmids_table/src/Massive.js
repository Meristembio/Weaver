import React from 'react'

class Massive extends React.Component {
    constructor(props) {
        super(props)
        this.state = {
            selected: [], // []
            plasmids_selected: [],
            mode: 'label', // ''
            newElementIDError: false,
            renderOutput: true
        }

        this.newElementID = React.createRef();
    }

    appendPlasmid(){

        let newValue = parseInt(this.newElementID.current.value)
        let newValidatedValue = false
        let newPlasmidElement = null

        this.props.plasmids.map((plasmid) => {
            if(parseInt(plasmid.ix) === newValue){
                newValidatedValue = newValue
                newPlasmidElement = plasmid
                newPlasmidElement.date = new Date().toISOString().slice(0, 10)
                newPlasmidElement.colony = null
                newPlasmidElement.conc = null
                return
            }
        })

        if (newValidatedValue){
            this.setState({
                selected: [...new Set([...this.state.selected, newValidatedValue])],
                plasmids_selected: [...this.state.plasmids_selected, newPlasmidElement],
                newElementIDError: false
            })
            this.newElementID.current.value = ''
        } else {
            this.setState({
                newElementIDError: true
            })
        }
    }

    removePlasmid(id) {
        this.setState({
            selected: this.state.selected.filter(item => item !== id)
        });
    }

    clearSelection() {
        this.setState({
            selected: []
        });
    }

    renderMassiveOutput() {
        if(this.state.mode === "label"){
            return <div id='labels'>
                    {this.state.plasmids_selected.map((plasmid, index) => {
                        return <div key={index} className="label-main">
                            <div className="label-info">
                                <p className="label-id-colony">
                                    <span class="label-value">ID {plasmid.ix} ~ c <span id="target-colony">{plasmid.colony ? plasmid.colony: "-"}</span></span>
                                </p>
                                <p className="label-name">
                                    {plasmid.n}
                                </p>
                                <p className="label-created">
                                    <span class="label-name">Date</span>
                                    <span class="label-value"><span className="target-date">{plasmid.date}</span></span>
                                </p>
                                {plasmid.conc && <p className="label-conc-quantus">
                                    <span class="label-name">Conc</span>
                                    <span class="label-value"><span className="target-concentration">{plasmid.conc}</span> [ng/ul]</span>
                                </p>}
                                <p className="label-size">
                                    <span class="label-name">Size</span>
                                    <span class="label-value">{plasmid.c} bp</span>
                                </p>
                            </div>
                        </div>
                    })}
            </div>
        }
    }

    renderMassiveInput() {
        return <div className='my-4'>
            <table id='plasmid-table-massive-selected' className="table sortable">
                <thead>
                    <tr>
                        <th>ID</th>
                        <th>Name</th>
                        {this.state.mode === "label" && [
                            <th key={1}>Date</th>,
                            <th key={2}>Colony</th>,
                            <th key={3}>Concentration</th>
                        ]}
                        <th>Remove</th>
                    </tr>
                </thead>
                <tbody>
                    {this.state.plasmids_selected.map((plasmid, index) => {
                        return <tr key={index}>
                        <td>
                            {plasmid.ix}
                            <input 
                                name={'id-'+plasmid.ix}
                                type="hidden" 
                                className="form-control" 
                                value={plasmid.ix}
                            />
                        </td>
                        <td>{plasmid.n}</td>
                        {this.state.mode === "label" && [
                            <td key={1}>
                                <input 
                                    name={'date-'+plasmid.ix}
                                    type="date" 
                                    className="form-control" 
                                    aria-label="Date"
                                    value={plasmid.date}
                                    onChange={(e)=>{
                                        this.setState(prev => ({
                                            plasmids_selected: prev.plasmids_selected.map(p =>
                                                String(p.ix) === String(plasmid.ix) ? { ...p, date: e.target.value } : p
                                            )
                                        }));
                                    }}
                                />
                            </td>,
                            <td key={2}>
                                <input 
                                    name={'colony-'+plasmid.ix}
                                    type="number" 
                                    className="form-control" 
                                    aria-label="Colony" 
                                    min="1"
                                    defaultValue={plasmid.wc}
                                    ref={plasmid.colony}
                                />
                            </td>,
                            <td key={3}>
                                <input 
                                    name={'conc-'+plasmid.ix}
                                    type="number" 
                                    className="form-control" 
                                    aria-label="Concentration" 
                                    min="1"
                                    ref={plasmid.conc}
                                />
                            </td>,
                        ]}
                        <td>
                            <button 
                                className='btn btn-outline-danger' 
                                onClick={() => this.removePlasmid(parseInt(plasmid.ix))}
                            >
                                <i className="bi bi-trash"></i>
                            </button>
                        </td>
                    </tr>
                    })}
                </tbody>
            </table>
            <button type="submit" className="btn btn-primary" onClick={()=>{
                this.setState({renderOutput:true})
            }}>Perform</button>
        </div>
    }

    render() {
        return (
            <div id='plasmid-massive'>
                <div className='row'>
                    <div className='col-md-3'>
                        <div id='plasmid-table-massive-select'>
                            <h4>Select by plasmid ID</h4>
                            <div className="input-group">
                                <input type="text" className={this.state.newElementIDError ? "form-control is-invalid": "form-control"} placeholder="Plasmid ID" ref={this.newElementID} onKeyDown={(e) => {
                                    if(e.key === "Enter") this.appendPlasmid()
                                }} />
                                <button className='btn btn-primary' data-bs-original-title="Append to selection" aria-label="Append plasmid" data-bs-toggle="tooltip" data-bs-placement="top" onClick={() => {this.appendPlasmid()}}><i className="bi bi-arrow-right"></i></button>
                                <button type="button" className="btn btn-outline-secondary ms-2" data-bs-original-title="Clear" aria-label="Clear" data-bs-toggle="tooltip" data-bs-placement="top" onClick={() => {this.clearSelection()}}><i className="bi bi-x-lg"></i></button>
                            </div>
                        </div>
                    </div>
                    <div className='col-md-9 d-flex flex-row-reverse'>
                        <div id="plasmids-filter" className="pe-table-filter">
                            <div className="pe-table-filter-buttons">    
                                <div className="pe-table-filter-buttons-group">
                                    <div className="pe-table-filter-buttons-group-title">
                                        Application
                                    </div>
                                    <div className="pe-table-filter-buttons-group-ops">                                        
                                        <div id='plasmid-table-massive-modes' className='mb-2'>
                                            <div className="btn-group" role="group" aria-label="Massive modes">
                                                <input type="radio" className="btn-check" name="massive_mode_radio" id="massive_mode-1" autoComplete="off" onChange={() => {this.setState({mode:'label'})}} />
                                                <label className="btn btn-outline-primary" htmlFor="massive_mode-1">Label</label>

                                                <input type="radio" className="btn-check" name="massive_mode_radio" id="massive_mode-2" autoComplete="off" onChange={() => {this.setState({mode:'validate'})}}/>
                                                <label className="btn btn-outline-primary" htmlFor="massive_mode-2">Validate</label>                                
                                            </div>
                                        </div>         
                                    </div>
                                </div>    
                            </div>
                        </div>
                    </div>
                </div>
                <div className='row'>
                    <div className='col'>
                        {this.state.mode === "" ? 
                        <div className='alert alert-info my-4'>Select an application to continue</div>:
                        this.state.selected.length > 0 ?  
                            this.state.renderOutput && [this.renderMassiveOutput(),
                                this.renderMassiveInput()]:
                                <div className='alert alert-info my-4'>Select plasmids to continue</div>
                        }
                    </div>
                </div>
            </div>
        );
    }
}

export default Massive;
