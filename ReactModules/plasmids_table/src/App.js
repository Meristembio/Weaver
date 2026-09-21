import './App.css';
import React from 'react'
//import Massive from './Massive';
import PlasmidsTable from './PlasmidsTable';
// comment on production mode
// import testData from './test.json'

const test_mode = false
// comment on test mode
const testData = []
const RECENTLY_VIEWED_STORAGE_KEY = 'weaver.recently-viewed-plasmids'
const MAX_RECENTLY_VIEWED = 20

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
        if (input) {
            input.addEventListener('input', this.handleSearchInput)
        }
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

    getRecentlyViewedIds = () => {
        try {
            const stored = window.localStorage.getItem(RECENTLY_VIEWED_STORAGE_KEY)
            const parsed = stored ? JSON.parse(stored) : []
            return Array.isArray(parsed) ? parsed.filter((id) => typeof id === 'string') : []
        } catch (error) {
            return []
        }
    }

    handlePlasmidView = (plasmidId) => {
        const id = String(plasmidId)
        const recentIds = [id, ...this.getRecentlyViewedIds().filter((storedId) => storedId !== id)]
            .slice(0, MAX_RECENTLY_VIEWED)
        try {
            window.localStorage.setItem(RECENTLY_VIEWED_STORAGE_KEY, JSON.stringify(recentIds))
        } catch (error) {
            // Browsers may disable localStorage; the plasmid link should still work.
        }
        this.loadPage(this.state.page, this.state.search, this.state.searchField)
    }

    clearRecentlyViewed = () => {
        try {
            window.localStorage.removeItem(RECENTLY_VIEWED_STORAGE_KEY)
        } catch (error) {
            // Browsers may disable localStorage; the list can still be refreshed.
        }
        this.loadPage(this.state.page, this.state.search, this.state.searchField)
    }

    loadPage = (page = this.state.page, search = this.state.search, searchField = this.state.searchField) => {
        const axios = require('axios')
        if (test_mode) {
            this.setState({ data: testData })
            return
        }
        const requestNumber = ++this.requestNumber
        this.setState({ loading: true })
        axios.get('/inventory/api/plasmids/', {
            params: {
                page: page,
                page_size: 50,
                q: search,
                search: searchField,
                recent_ids: this.getRecentlyViewedIds().join(',')
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

    renderRecentlyViewed = () => {
        const recentlyViewed = this.state.data.recently_viewed || []
        if (!recentlyViewed.length) return null
        const recentlyViewedData = { ...this.state.data, plasmids: recentlyViewed }
        return <section className="mb-4" aria-labelledby="recently-viewed-heading">
            <div className="d-flex align-items-center mb-3">
                <h4 id="recently-viewed-heading" className="fw-bold mb-0">Recently viewed</h4>
                <button type="button" className="btn btn-outline-secondary btn-sm ms-2 p-1"
                    onClick={this.clearRecentlyViewed}
                    aria-label="Clear recently viewed">
                    <i className="bi bi-trash3"></i>
                </button>
            </div>
            <PlasmidsTable
                data={recentlyViewedData}
                tableId="recently-viewed-plasmids-table"
                tablePrefix="recently-viewed"
                onPlasmidView={this.handlePlasmidView}
            />
        </section>
    }

    renderPaginationControls = () => {
        const pageCount = this.state.data.num_pages || 1
        const firstPageDisabled = this.state.page <= 1 || this.state.loading
        const lastPageDisabled = this.state.page >= pageCount || this.state.loading
        return <div className="d-flex justify-content-between align-items-center my-2">
            <small className="text-secondary">
                Showing {this.state.data.plasmids.length} of {this.state.data.total} plasmids
            </small>
            <nav aria-label="Plasmid pagination">
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
        if (this.state.data) {
            if (this.state.data.plasmids) {
                const recentlyViewed = this.renderRecentlyViewed()
                const plasmidList = <>
                    <h4 className="fw-bold mb-3">Plasmid List</h4>
                    {this.renderPaginationControls()}
                    <PlasmidsTable data={this.state.data} onPlasmidView={this.handlePlasmidView} />
                    {this.renderPaginationControls()}
                </>
                return <div id='plasmid-wrapper'>
                    <div id='plasmid-massive-wrapper'>
                        {/*<Massive plasmids={this.state.data.plasmids} />*/}
                    </div>
                    <div id='plasmid-table-wrapper'>
                        {this.state.search ? <>{plasmidList}{recentlyViewed}</> : <>{recentlyViewed}{plasmidList}</>}
                    </div>
                </div>
            } else {
                return <div className="alert alert-info">
                    <i className="bi bi-emoji-frown"></i> No plasmids
                </div>
            }
        } else {
            return <div className="alert alert-info">
                <div className="spinner-grow spinner-grow-sm" role="status">
                    <span className="visually-hidden">...</span>
                </div> Loading plasmids
            </div>
        }
    }
}

export default App;
