import './App.css';
import React from 'react'

function StockRow({ glycerolstock }) {
    const plasmidName = glycerolstock.pn || ''
    const plasmidIndex = glycerolstock.pix || ''
    const stockName = plasmidName
        ? `${glycerolstock.s} / ${plasmidName}`
        : glycerolstock.s
    const searchAll = [
        glycerolstock.s,
        plasmidName,
        plasmidIndex,
        glycerolstock.bn,
        glycerolstock.bl,
        glycerolstock.qr_id
    ].filter(Boolean).join(' ')
    const searchName = [glycerolstock.s, plasmidName].filter(Boolean).join(' ')

    let plasmidOutput = <span>No plasmid</span>
    if (glycerolstock.pi) {
        plasmidOutput = <a href={`/inventory/plasmid/${glycerolstock.pi}`} className="btn btn-outline-secondary"
            role="button">
            <span>{plasmidName}</span>
            <span className="plasmid_list-id badge text-bg-light border ms-1">{plasmidIndex}</span>
        </a>
    }

    let editOutput = null
    if (glycerolstock.p) {
        editOutput = <a href={`/inventory/glycerolstock/edit/${glycerolstock.i}`}
            className="btn text-secondary" role="button" data-bs-toggle="tooltip"
            data-bs-placement="top" title="Edit"><i className="bi bi-pencil-fill"></i></a>
    }

    let plasmidIcon = null
    if (glycerolstock.pcs === 'v') {
        plasmidIcon = <i className="bi bi-check-circle ms-2" data-bs-toggle="tooltip"
            data-bs-placement="top" title="Validated"></i>
    } else if (glycerolstock.pcs === 'r') {
        plasmidIcon = <i className="bi bi-bookmarks ms-2" data-bs-toggle="tooltip"
            data-bs-placement="top" title="Reference"></i>
    } else if (glycerolstock.pcs === 'c') {
        plasmidIcon = <i className="bi bi-hammer ms-2" data-bs-toggle="tooltip"
            data-bs-placement="top" title="Under construction"></i>
    }

    let filterClasses = 'filter-item'
    if (glycerolstock.pl != null && glycerolstock.pl !== '') {
        filterClasses += ` filter-l${glycerolstock.pl}`
    }
    if (glycerolstock.pt != null && glycerolstock.pt !== '') {
        filterClasses += ` filter-t${glycerolstock.pt}`
    }

    return <tr className={filterClasses}>
        <td>
            <a className="btn btn-success table-search-search_on me-1"
                data-search-all={searchAll}
                data-search-name={searchName}
                data-search-idx={plasmidIndex}
                data-name={stockName}
                role="button" href={`/inventory/glycerolstock/${glycerolstock.i}`}>
                <span>{stockName}</span>
                {plasmidIndex && <span className="plasmid_list-id badge text-bg-light text-success ms-1">{plasmidIndex}</span>}
            </a>
        </td>
        <td>
            {editOutput}
            <a href={`/inventory/glycerolstock/label/${glycerolstock.i}`} className="btn text-info me-1"
                role="button" data-bs-toggle="tooltip" data-bs-placement="top" title="Print label">
                <i className="bi bi-tag-fill"></i>
            </a>
        </td>
        <td>{glycerolstock.s}</td>
        <td>{plasmidOutput}{plasmidIcon}</td>
        <td>{glycerolstock.br}{glycerolstock.bc}</td>
        <td>{glycerolstock.bn}</td>
        <td>{glycerolstock.bl}</td>
    </tr>
}

function StocksTable({ data }) {
    return <table id="glycerolstocks-table" className="table table-striped table-hover table-search-target align-middle">
        <thead>
            <tr>
                <th scope="col">Stock</th>
                <th scope="col">Actions</th>
                <th scope="col">Strain</th>
                <th scope="col">Plasmid</th>
                <th scope="col">Position</th>
                <th scope="col">Box</th>
                <th scope="col">Location</th>
            </tr>
        </thead>
        <tbody>
            {data.glycerolstocks.map((glycerolstock) => <StockRow
                key={glycerolstock.i} glycerolstock={glycerolstock} />)}
        </tbody>
    </table>
}

class App extends React.Component {
    constructor(props) {
        super(props)
        this.state = {
            data: false,
            page: 1,
            search: '',
            searchField: 'all',
            loading: false
        }
        this.searchTimeout = null
        this.requestNumber = 0
    }

    componentDidMount() {
        this.loadPage()
        const input = document.getElementById('table_search-input')
        if (input) input.addEventListener('input', this.handleSearchInput)
        document.querySelectorAll('.table_search-target').forEach((button) => {
            button.addEventListener('click', this.handleSearchTarget)
        })
        const clearButton = document.getElementById('table_search-clear')
        if (clearButton) clearButton.addEventListener('click', this.handleSearchClear)
    }

    componentWillUnmount() {
        const input = document.getElementById('table_search-input')
        if (input) input.removeEventListener('input', this.handleSearchInput)
        document.querySelectorAll('.table_search-target').forEach((button) => {
            button.removeEventListener('click', this.handleSearchTarget)
        })
        const clearButton = document.getElementById('table_search-clear')
        if (clearButton) clearButton.removeEventListener('click', this.handleSearchClear)
        clearTimeout(this.searchTimeout)
    }

    componentDidUpdate() {
        window.onReady()
    }

    loadPage = (page = this.state.page, search = this.state.search, searchField = this.state.searchField) => {
        const axios = require('axios')
        const requestNumber = ++this.requestNumber
        this.setState({ loading: true })
        axios.get('/inventory/api/glycerolstocks/', {
            params: {
                page: page,
                page_size: 50,
                q: search,
                search: searchField
            }
        }).then((response) => {
            if (requestNumber !== this.requestNumber || !response.data) return
            this.setState({ data: response.data, page: response.data.page, loading: false })
        }).catch(() => {
            if (requestNumber === this.requestNumber) this.setState({ loading: false })
        })
    }

    handleSearchInput = (event) => {
        const search = event.target.value
        clearTimeout(this.searchTimeout)
        this.searchTimeout = setTimeout(() => {
            this.setState({ search, page: 1 })
            this.loadPage(1, search, this.state.searchField)
        }, 250)
    }

    handleSearchTarget = (event) => {
        const searchField = event.currentTarget.dataset.type || 'all'
        this.setState({ searchField, page: 1 })
        this.loadPage(1, this.state.search, searchField)
    }

    handleSearchClear = () => {
        const input = document.getElementById('table_search-input')
        if (input) input.value = ''
        clearTimeout(this.searchTimeout)
        this.setState({ search: '', page: 1 })
        this.loadPage(1, '', this.state.searchField)
    }

    changePage = (page) => {
        if (page < 1 || page > (this.state.data.num_pages || 1)) return
        this.setState({ page })
        this.loadPage(page, this.state.search, this.state.searchField)
    }

    renderPaginationControls = () => {
        const pageCount = this.state.data.num_pages || 1
        const firstPageDisabled = this.state.page <= 1 || this.state.loading
        const lastPageDisabled = this.state.page >= pageCount || this.state.loading
        return <div className="d-flex justify-content-between align-items-center my-2">
            <small className="text-secondary">
                Showing {this.state.data.glycerolstocks.length} of {this.state.data.total} stocks
            </small>
            <nav aria-label="Stock pagination">
                <ul className="pagination justify-content-end mb-0 plasmid-pagination">
                    <li className={'page-item' + (firstPageDisabled ? ' disabled' : '')}>
                        <button className="page-link" disabled={firstPageDisabled}
                            onClick={() => this.changePage(1)} aria-label="First page">&lt;&lt;</button>
                    </li>
                    <li className={'page-item' + (firstPageDisabled ? ' disabled' : '')}>
                        <button className="page-link" disabled={firstPageDisabled}
                            onClick={() => this.changePage(this.state.page - 1)} aria-label="Previous page">Previous</button>
                    </li>
                    <li className="page-item disabled"><span className="page-link">{this.state.page} / {pageCount}</span></li>
                    <li className={'page-item' + (lastPageDisabled ? ' disabled' : '')}>
                        <button className="page-link" disabled={lastPageDisabled}
                            onClick={() => this.changePage(this.state.page + 1)} aria-label="Next page">Next</button>
                    </li>
                    <li className={'page-item' + (lastPageDisabled ? ' disabled' : '')}>
                        <button className="page-link" disabled={lastPageDisabled}
                            onClick={() => this.changePage(pageCount)} aria-label="Last page">&gt;&gt;</button>
                    </li>
                </ul>
            </nav>
        </div>
    }

    render() {
        if (!this.state.data) {
            return <div className="alert alert-info">
                <div className="spinner-grow spinner-grow-sm" role="status">
                    <span className="visually-hidden">...</span>
                </div> Loading glycerol stocks
            </div>
        }

        if (!this.state.data.total) {
            return <div className="alert alert-warning">No glycerol stocks found.</div>
        }

        return <>
            {this.renderPaginationControls()}
            <StocksTable data={this.state.data} />
            {this.renderPaginationControls()}
        </>
    }
}

export default App
