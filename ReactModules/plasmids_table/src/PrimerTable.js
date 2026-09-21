import React from 'react';

function tmValue(sequence) {
    const counts = { a: 0, c: 0, g: 0, t: 0 };
    sequence.toLowerCase().split('').forEach((base) => {
        if (Object.prototype.hasOwnProperty.call(counts, base)) counts[base] += 1
    })
    const gc = counts.g + counts.c
    const at = counts.a + counts.t
    if (sequence.length < 14) return (at * 2) + (gc * 4)
    return 64.9 + (41 * (gc - 16.4) / (at + gc))
}

function PrimerRow({ primer }) {
    const completeSequence = (primer.sequence_5 || '') + (primer.sequence_3 || '')
    const primerId = primer.display_idx || ''
    const sequence3Length = (primer.sequence_3 || '').length
    const completeTm = tmValue(completeSequence).toFixed(1)
    const sequence3Tm = tmValue(primer.sequence_3 || '').toFixed(1)
    return <tr>
        <td>
            <a className="btn-group table-search-search_on me-1" role="button"
                href={`/inventory/primer/${primer.id}/`}
                data-search-all={`${primerId} ${primer.display_name} ${primer.intended_use} ${completeSequence}`}
                data-search-idx={primerId}
                data-search-name={primer.display_name}>
                {primerId && <span className="btn btn-success opacity-75 fw-bold"
                    data-bs-toggle="tooltip" data-bs-placement="top"
                    data-bs-custom-class="custom-tooltip" data-bs-title="Primer ID">{primerId}</span>}
                <span className="btn btn-success">{primer.display_name}</span>
            </a>
            <a href={`/inventory/primer/edit/${primer.id}/`}
                className={'btn text-secondary' + (primer.can_edit ? '' : ' disabled')}
                role="button"><i className="bi bi-pencil-fill"></i></a>
        </td>
        <td><span className="overhang">{primer.sequence_5}</span>{primer.sequence_3}</td>
        <td>{completeSequence.length} / {sequence3Length}</td>
        <td>{completeTm} / {sequence3Tm} °C</td>
    </tr>
}

function PrimerTable({ data }) {
    return <table className="table table-striped table-hover table-search-target align-middle">
        <thead>
            <tr>
                <th scope="col">Primer</th>
                <th scope="col">Sequence (5' → 3') (5' overhang underlined)</th>
                <th scope="col">Length (total / 3')</th>
                <th scope="col">Tm (total / 3')</th>
            </tr>
        </thead>
        <tbody>
            {data.primers.map((primer) => <PrimerRow key={primer.id} primer={primer} />)}
        </tbody>
    </table>
}

export default PrimerTable;
