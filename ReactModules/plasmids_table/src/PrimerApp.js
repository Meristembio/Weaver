import './App.css';
import React from 'react';
import PrimerTable from './PrimerTable';

class PrimerApp extends React.Component {
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
        axios.get('/inventory/api/primers/', {
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
                Showing {this.state.data.primers.length} of {this.state.data.total} primers
            </small>
            <nav aria-label="Primer pagination">
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
                </div> Loading primers
            </div>
        }

        if (!this.state.data.total) {
            return <div className="alert alert-warning">No primers found.</div>
        }

        return <>
            {this.renderPaginationControls()}
            <PrimerTable data={this.state.data} />
            {this.renderPaginationControls()}
        </>
    }
}

export default PrimerApp;
