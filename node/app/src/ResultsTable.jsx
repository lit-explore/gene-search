import React, { useMemo } from 'react';
import PropTypes from 'prop-types';
import BTable from 'react-bootstrap/Table';
import { useTable, usePagination } from 'react-table'

const ResultsTable = (props) => {
  // table columns
  const columns = useMemo(() => [
    {
    // https://pubmed.ncbi.nlm.nih.gov/24949858/
      Header: 'PMID',
      accessor: 'id',
      Cell: function addLink({ row }) {
        const url = `https://pubmed.ncbi.nlm.nih.gov/${row.original.id}/`
        const html = `<a target='blank' href='${url}'>${row.original.id}</a>`

        return(
          <span dangerouslySetInnerHTML={{ __html: html }} /> 
        )
      }
    },
    {
      Header: 'Citation',
      accessor: 'citation',
    },
    {
      Header: 'Genes',
      accessor: 'genes',
    },
    {
      Header: '# Genes',
      accessor: 'num_genes',
      width: 250
    },
    {
      Header: 'Cluster',
      accessor: 'cluster',
    }
  ], [])

  // generate table data array
  const data = useMemo(() => props.data, [props.data])

  const { 
    getTableProps,
    getTableBodyProps,
    headerGroups,
    prepareRow,
    pageOptions,
    page,
    state: { pageIndex, pageSize },
    gotoPage,
    previousPage,
    nextPage,
    setPageSize,
    canPreviousPage,
    canNextPage,
  } = useTable({
      columns,
      data,
      initialState: { pageIndex: 2, pageSize: 5 },
   },
    usePagination
  )

  const pageSizeOptions = [5].concat([10, 25, 50, 100].filter(x => x <= props.data.length))

  // memoize table html to avoid re-rendering on hover state changes..
  if (props.data.length === 0) {
    return(<div/>)
  }

  return (
      <div id='results-table'>
        <BTable striped bordered hover size="sm" {...getTableProps()}>
          <thead>
            {headerGroups.map((headerGroup, headerInd) => (
              <tr key={"header-" + headerInd} {...headerGroup.getHeaderGroupProps()}>
                {headerGroup.headers.map((column, colInd) => (
                  <th key={"col-" + colInd} {...column.getHeaderProps()}>
                    {column.render('Header')}
                  </th>
                ))}
              </tr>
            ))}
          </thead>
          <tbody {...getTableBodyProps()}>
            {page.map((row, i) => {
              prepareRow(row)
              return (
                <tr style={{'color': row.original.color}} id={"table-row-" + row.original.id} key={"row-" + i} {...row.getRowProps()}>
                  {row.cells.map((cell, j) => {
                    return (
                      <td key={"cell-" + j} {...cell.getCellProps()}>
                        {cell.render('Cell')}
                      </td>
                    )
                  })}
                </tr>
              )
            })}
          </tbody>
        </BTable>
        <div>
          <button onClick={() => previousPage()} disabled={!canPreviousPage}>
            Prev
          </button>
          <button onClick={() => nextPage()} disabled={!canNextPage}>
            Next
          </button>
          <span>
            Page{' '}
            <em>
              {pageIndex + 1} of {pageOptions.length}
            </em>
          </span> || 
          <span>Go to: </span>
          <input
            type="number"
            style={{"width": "75px"}}
            defaultValue={pageIndex + 1 || 1}
            onChange={e => {
              const page = e.target.value ? Number(e.target.value) - 1 : 0
              gotoPage(page)
            }}
          />
          <select
            value={pageSize}
            onChange={e => {
              setPageSize(Number(e.target.value))
            }}
          >
            {pageSizeOptions.map(pageSize => (
              <option key={pageSize} value={pageSize}>
                Show {pageSize}
              </option>
            ))}
          </select>
        </div>
      </div>
    )
}

ResultsTable.propTypes = {
  data: PropTypes.array,
};

export default ResultsTable;
