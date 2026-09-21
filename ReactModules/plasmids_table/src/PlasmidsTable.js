import React from 'react'

class PlasmidElement extends React.Component {
    render() {
        const plasmid = this.props.plasmid
        const tablePrefix = this.props.tablePrefix || 'plasmids'
        const elementPrefix = tablePrefix + '-' + plasmid.i
        let output = []

        let plasmid_level = ""
        if (plasmid.l !== null) plasmid_level = " filter-l" + plasmid.l
        let plasmid_type = ""
        if (plasmid.t !== null) plasmid_type = " filter-t" + plasmid.t
        let plasmid_part = ""
        if (plasmid.pk) plasmid_part = " filter-ap-" + plasmid.pk.replaceAll("_", "-")

        let plasmid_computed_size = 'No edit perms'
        let plasmid_create_build = 'No edit perms'

        if (plasmid.p) {
            let reOpts = []
            this.props.RESTRICTION_ENZYMES.forEach((re) => {
                let btnClass = "btn-outline-primary"
                if (plasmid.r.toLowerCase() === re.name.toLowerCase()) {
                    btnClass = "btn-primary"
                }
                reOpts.push(
                    <button
                        key={re.name}
                        className={"btn " + btnClass + " btn-sm me-1"}
                        role="button" name="enzyme" value={re.name}>
                        {re.name}
                    </button>
                )
            })
            plasmid_create_build = <div>
                <form method="post" className="default-style inline" target="_blank"
                    action={"/inventory/plasmid/view_edit/" + plasmid.i}>
                    <input type="hidden" name="csrfmiddlewaretoken" value={this.props.csrf_token} />
                    <button className="btn text-success me-1" name="create" data-bs-toggle="tooltip" data-bs-placement="top" title="Create blank"><i
                        className="bi bi-file-earmark-plus"></i></button>
                </form>
                <div className="dropdown dropdown-enzymes" data-bs-toggle="tooltip" data-bs-placement="top" title="Build">
                    <button type="button" className="btn text-primary dropdown-toggle" id={elementPrefix + "-dropdownEnzymes"} data-bs-toggle="dropdown" aria-expanded="false"><i
                        className="bi bi-hammer"></i></button>
                    <div className="dropdown-menu p-2 fw-light " aria-labelledby={elementPrefix + "-dropdownEnzymes"}>
                        <div className="dropdown-menu-header">Chooose enzyme</div>
                        <hr className="m-1" />
                        <form method="post" className="default-style" target="_blank">
                            <input type="hidden" name="csrfmiddlewaretoken" value={this.props.csrf_token} />
                            <input type="hidden" name="create_from_parts" />
                            <div className="dropdown-menu-body pt-1">
                                {reOpts}
                            </div>
                        </form>
                    </div>
                </div>
            </div>
        }

        let plasmid_insert_computed_size = ""

        let plasmid_edits = []
        if (plasmid.p) {
            plasmid_edits.push(<a
                key={plasmid.i}
                href={"/inventory/plasmid/validation/edit/" + plasmid.i}
                className="btn text-primary me-1"
                role="button" target="_blank" rel="noreferrer" data-bs-toggle="tooltip" data-bs-placement="top" title="Edit validation"><i
                    className="bi bi-check2-square"></i></a>)
        }
        let plasmid_sequence_options = []

        if (plasmid.hs) {
            plasmid_sequence_options.push(<a key={plasmid.i + "ve"} href={"/inventory/plasmid/view_edit/" + plasmid.i} className="btn text-info me-1"
                role="button" target="_blank" rel="noreferrer" data-bs-toggle="tooltip" data-bs-placement="top" title="View / edit sequence"><i className="bi bi-eye"></i> / <i
                    className="bi bi-pencil"></i></a>)
            plasmid_sequence_options.push(
                <div key={plasmid.i + "dl"} className="dropdown dropdown-download" data-bs-toggle="tooltip"
                    data-bs-placement="top" title="Download">
                        <button type="button" className="btn text-primary dropdown-toggle me-1"
                        id={elementPrefix + "-dropdownDownload"} data-bs-toggle="dropdown" aria-expanded="false"><i
                            className="bi bi-download"></i>
                    </button>
                    <div className="dropdown-menu p-2 fw-light " aria-labelledby={elementPrefix + "-dropdownDownload"}>
                        <div className="dropdown-menu-header">Chooose format</div>
                        <hr className="m-1" />
                        <div className="dropdown-menu-body pt-1">
                            <a href={"/inventory/plasmid/download/" + plasmid.i}
                                className="btn btn-outline-primary btn-sm me-1"
                                role="button"
                                download={plasmid.n}>ORIG</a>
                            <a href={"/inventory/plasmid/download/" + plasmid.i + "?format=gb"}
                                className="btn btn-outline-primary btn-sm me-1"
                                role="button"
                                download={plasmid.n}>GB</a>
                            <a href={"/inventory/plasmid/download/" + plasmid.i + "?format=fasta"}
                                className="btn btn-outline-primary btn-sm" role="button"
                                download={plasmid.n}>FASTA</a>
                        </div>
                    </div>
                </div>
            )
            plasmid_sequence_options.push(<a key={plasmid.i + "tl"} href={"/inventory/plasmid/digest/" + plasmid.i}
                className="btn text-secondary me-1" data-bs-toggle="tooltip" data-bs-placement="top" title="Digest"
                role="button" target="_blank" rel="noreferrer"><i className="bi bi-scissors"></i></a>)
            plasmid_sequence_options.push(<a key={plasmid.i + "dp"} href={"/inventory/plasmid/pcr/" + plasmid.i} className="btn text-success me-1"
                role="button" target="_blank" rel="noreferrer" data-bs-toggle="tooltip" data-bs-placement="top" title="Design PCR"><i
                    className="bi bi-arrow-return-right"></i></a>)
            plasmid_sequence_options.push(<div key={plasmid.i + "da"} className="dropdown dropdown-align" data-bs-toggle="tooltip" data-bs-placement="top" title="Align">
                <button type="button" className="btn text-warning dropdown-toggle" id={elementPrefix + "-dropdownAlign"} data-bs-toggle="dropdown" aria-expanded="false"><i className="bi bi-list-nested"></i>
                </button>
                <div className="dropdown-menu p-2 fw-light " aria-labelledby={elementPrefix + "-dropdownAlign"}>
                    <div className="dropdown-menu-header">Chooose type</div>
                    <hr className="m-1" />
                    <div className="dropdown-menu-body pt-1">
                        <a href={"/inventory/plasmid/align/fasta/" + plasmid.i} className="btn btn-outline-primary btn-sm me-1" role="button" target="_blank">Fasta</a>
                        <a href={"/inventory/plasmid/align/sanger/" + plasmid.i} className="btn btn-outline-primary btn-sm me-1" role="button" target="_blank">Sanger</a>
                    </div>
                </div>
            </div>)

            if (plasmid.c !== null && plasmid.c > 0) {
                if (plasmid.c < 1000) {
                    plasmid_computed_size = plasmid.c + " b"
                } else {
                    plasmid_computed_size = Math.floor(plasmid.c / 1000 * 10) / 10 + " k"
                }

            }
            let plasmid_size_output = [plasmid_computed_size]
            if (plasmid.ic !== null && plasmid.ic > 0) {
                if (plasmid.ic < 1000) {
                    plasmid_insert_computed_size = plasmid.ic + " b"
                } else {
                    plasmid_insert_computed_size = Math.floor(plasmid.ic / 1000 * 10) / 10 + " k"
                }
                plasmid_size_output.push(<span key={plasmid.i}> / {plasmid_insert_computed_size}</span>)
            }
            plasmid_computed_size = <span>{plasmid_size_output}</span>
            plasmid_create_build = null
        }
        else {
            plasmid_computed_size = <span>No sequence<br />{plasmid_create_build}</span>
        }
        let table_filters_output = ""
        this.props.table_filters.forEach((table_filter) => {
            if (table_filter[0] === 'startswith') {
                table_filter[2].forEach((table_filter_op) => {
                    if (plasmid.n.toLowerCase().startsWith(table_filter_op[1])) {
                        table_filters_output += " filter-" + table_filter_op[1]
                    }
                })
            }
        })
        let plasmid_icon = ""
        if (plasmid.cs === 'v') {
            // verified
            plasmid_icon = <i className="bi bi-check-circle ms-2" data-bs-toggle="tooltip" data-bs-placement="top" title="Validated"></i>
        } else if (plasmid.cs === 'r') {
            // reference
            plasmid_icon = <i className="bi bi-bookmarks ms-2" data-bs-toggle="tooltip" data-bs-placement="top" title="Reference"></i>
        } else if (plasmid.cs === 'c') {
            // under construction
            plasmid_icon = <i className="bi bi-hammer ms-2" data-bs-toggle="tooltip" data-bs-placement="top" title="Under construction"></i>
        }
        let plasmid_type_output = []
        let plasmid_level_output = []
        let plasmid_part_output = []
        if (plasmid.l !== undefined) {
            plasmid_level_output.push(plasmid.l)
        }
        if (plasmid.t !== undefined) {
            let type_name = "Insert"
            if (plasmid.t === 1) {
                type_name = "Receiver"
            }
            plasmid_type_output.push(type_name)
        }
        if (plasmid.pnm) {
            plasmid_part_output.push(plasmid.pnm)
        }
        let plasmid_name = plasmid.n
        if (plasmid.n.length > 25)
            plasmid_name = plasmid.n.substring(0, 26) + "..."
        let plasmid_edit = ""
        if (plasmid.p)
            plasmid_edit = <a href={"/inventory/plasmid/edit/" + plasmid.i} className="btn text-secondary me-1"
                role="button" target="_blank" rel="noreferrer" data-bs-toggle="tooltip" data-bs-placement="top" title="Edit"><i className="bi bi-pencil-fill"></i></a>
        output.push(<tr key={plasmid.i} className={"filter-item" + plasmid_level + plasmid_type + plasmid_part + table_filters_output}>
            <td>
                <a href={"/inventory/plasmid/" + plasmid.i}
                    className="btn-group table-search-search_on me-1"
                    data-name={plasmid.n}
                    data-plasmid-id={plasmid.i}
                    data-search-all={plasmid.n + plasmid.ix + plasmid.d + plasmid.iu}
                    data-search-idx={plasmid.ix}
                    data-search-name={plasmid.n}
                    role="button" target="_blank" rel="noreferrer"
                    onClick={() => this.props.onPlasmidView(plasmid.i)}>
                    <button className="btn btn-success plasmid_list-name">{plasmid_name}</button>
                    <button className="btn btn-success opacity-75 plasmid_list-id fw-bold">{plasmid.ix}</button>
                </a>
                {plasmid_icon}
            </td>
            <td>
                {plasmid_edit}
                <button className="btn text-secondary me-1 copy_clipboard-child" data-bs-toggle="tooltip" data-bs-placement="top" title="Copy name + id"><i className="bi bi-clipboard copy_clipboard" data-cc={plasmid.cn}></i></button>
                <div className="dropdown dropdown-align" data-bs-toggle="tooltip" data-bs-placement="top" title="Copy options">
                    <button type="button" className="btn text-secondary dropdown-toggle" id={elementPrefix + "-dropdownCopy"} data-bs-toggle="dropdown" aria-expanded="false"><i className="bi bi-clipboard2-minus"></i>
                    </button>
                    <div className="dropdown-menu p-2 fw-light " aria-labelledby={elementPrefix + "-dropdownCopy"}>
                        <div className="dropdown-menu-header">Options</div>
                        <hr className="m-1" />
                        <div className="dropdown-menu-body pt-1">
                            <button className="dropdown-item btn text-secondary me-1 copy_clipboard-child" data-bs-toggle="tooltip" data-bs-placement="top" title="Copy name"><i className="bi bi-clipboard copy_clipboard" data-cc={plasmid.n}></i> Name</button>
                            <button className="dropdown-item btn text-secondary me-1 copy_clipboard-child" data-bs-toggle="tooltip" data-bs-placement="top" title="Copy name"><i className="bi bi-clipboard copy_clipboard" data-cc={plasmid.cn}></i> Name + ID</button>
                            <button className="dropdown-item btn text-secondary me-1 copy_clipboard-child" data-bs-toggle="tooltip" data-bs-placement="top" title="Copy name"><i className="bi bi-clipboard copy_clipboard" data-cc={plasmid.cn + " - Colony: " + plasmid.wc}></i> Name + ID + Colony</button>
                        </div>
                    </div>
                </div>
                <a href={"/inventory/plasmid/label/" + plasmid.i} className="btn text-info me-1"
                    role="button" target="_blank" rel="noreferrer" data-bs-toggle="tooltip" data-bs-placement="top" title="Print label"><i className="bi bi-tag-fill"></i></a>
                {plasmid_edits}
            </td>
            <td>{plasmid_create_build}{plasmid_sequence_options}</td>
            <td>
                <button type="button" className="btn btn-outline-secondary" data-bs-toggle="modal" data-bs-target={"#modal-" + tablePrefix + "-" + this.props.index}>
                    Details
                </button>
                <div className="modal" id={"modal-" + tablePrefix + "-" + this.props.index} tabIndex="-1" aria-labelledby={"modal_title-" + tablePrefix + "-" + this.props.index} aria-hidden="true">
                    <div className="modal-dialog">
                        <div className="modal-content">
                            <div className="modal-header">
                                <h1 className="modal-title fs-5" id={"modal_title-" + tablePrefix + "-" + this.props.index}>Details</h1>
                                <button type="button" className="btn-close" data-bs-dismiss="modal" aria-label="Close"></button>
                            </div>
                            <div className="modal-body">
                                <table className="table">
                                    <thead>
                                        <tr>
                                            <th scope="col">Property</th>
                                            <th scope="col">Value</th>
                                        </tr>
                                    </thead>
                                    <tbody>
                                        <tr>
                                            <td>Name</td>
                                            <td>{plasmid.n}</td>
                                        </tr>
                                        <tr>
                                            <td>ID</td>
                                            <td>{plasmid.ix}</td>
                                        </tr>
                                        <tr>
                                            <td>Level</td>
                                            <td>{plasmid_level_output}</td>
                                        </tr>
                                        <tr>
                                            <td>Type</td>
                                            <td>{plasmid_type_output}</td>
                                        </tr>
                                        <tr>
                                            <td>Part</td>
                                            <td>{plasmid_part_output}</td>
                                        </tr>
                                        <tr>
                                            <td>Marker</td>
                                            <td>{plasmid.sm}</td>
                                        </tr>
                                        <tr>
                                            <td>Total / insert length</td>
                                            <td>{plasmid_computed_size}</td>
                                        </tr>
                                        <tr>
                                            <td>Colony</td>
                                            <td>{plasmid.wc}</td>
                                        </tr>
                                        <tr>
                                            <td>Ligation concentration</td>
                                            <td>{plasmid.lc}</td>
                                        </tr>
                                        <tr>
                                            <td>Intended use</td>
                                            <td>{plasmid.iu}</td>
                                        </tr>
                                        <tr>
                                            <td>Description</td>
                                            <td>{plasmid.d}</td>
                                        </tr>
                                    </tbody>
                                </table>
                            </div>
                        </div>
                    </div>
                </div>
            </td>
        </tr>)
        return output
    }
}

class PlasmidsTable extends React.Component {
    render() {
        let output = []
        this.props.data.plasmids.forEach((plasmid, index) => {
            output.push(
                <PlasmidElement
                    key={index}
                    index={index}
                    plasmid={plasmid}
                    tablePrefix={this.props.tablePrefix || this.props.tableId || 'plasmids'}
                    onPlasmidView={this.props.onPlasmidView}
                    csrf_token={this.props.data.csrf_token}
                    table_filters={this.props.data.table_filters}
                    RESTRICTION_ENZYMES={this.props.data.RESTRICTION_ENZYMES}
                />)
        })
        return <table id={this.props.tableId || "plasmids-table"} className="table table-striped table-hover sortable table-search-target align-middle">
            <thead>
                <tr>
                    <th scope="col" data-defaultsort='disabled'>Plasmid</th>
                    <th scope="col" data-defaultsort='disabled'>Actions</th>
                    <th scope="col" data-defaultsort='disabled'>Sequence</th>
                    <th scope="col" data-defaultsort='disabled'>Details</th>
                </tr>
            </thead>
            <tbody>{output}</tbody>
        </table>
    }
}

export default PlasmidsTable;
