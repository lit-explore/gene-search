//
// Pubmed x Gene Search
// KH (Oct21)
//
import { useEffect, useState } from "react";
import { Button, Col, Container, Row, Form  } from "react-bootstrap";
import axios from "axios";
import ResultsTable from "./ResultsTable";
import Network from "./Network";
import useWindowSize from './hooks/useWindowSize';

function App() { 

  const [articleMatches, setArticleMatches] = useState([])
  const [networkData, setNetworkData] = useState({nodes: [], links: []})
  const [networkDims, setNetworkDims] = useState({width: 640, height: 480})

  const windowSize = useWindowSize();

  // set initial network dimensions, when windowSize hook is ready..
  useEffect(() => {
    if ((windowSize.height === undefined) || (windowSize.width === undefined)) {
      return
    }
    console.log("Window size changed:")
    console.log(windowSize.width)
    console.log(windowSize.height)

    setNetworkDims({
      width: windowSize.width * 0.65,
      height: windowSize.height * 0.5
    })
  }, [windowSize.height, windowSize.width])

  function handleSubmit(event) {
    // construct query URL
    const input = document.getElementById("gene-input").value 
    const genes = input.trim().replace(/[,\n]/g, ",")

    const limit = document.getElementById("article-limit").value

    const baseURL = "http://localhost:5000/query/"
    const queryURL = `${baseURL}?genes=${genes}&limit=${limit}`

    console.log(queryURL)

    // submit query to API
    axios.get(queryURL).then(resp => {
      console.log(resp)
      setArticleMatches(resp.data.articles)

      let articles = resp.data.articles

      // network tooltip text (title + genes)
      articles.forEach(entry => {
        entry.label = entry.citation + "<br /><b>" + entry.genes + "</b>"
      })

      // build graph data object
      setNetworkData({
        "nodes": articles,
        "links": resp.data.network.map(x => ({'source': x[0], 'target': x[1]}))
      })
    }).catch(function (error) {
      console.log("[Error] Error encountered!");
      console.log(error);
    });

    event.preventDefault()
  }

  return (
    <Container fluid>
      <Row>
    <h2 style={{margin: "20px 0px"}}>Pubmed Gene Set Search</h2>
      </Row>
      <Row>
        <Col xs={3}> 
          <Form onSubmit={handleSubmit}>
            <Form.Label>Gene Symbols (comma- or newline-separated)</Form.Label>
            <Form.Control id="gene-input" as="textarea" placeholder="Gene Symbols" style={{ height: '100px' }} />
            <Form.Select id="article-limit" aria-label="Article Limit">
              <option value="25">25</option>
              <option value="50">50</option>
              <option value="100">100</option>
              <option value="250">250</option>
              <option value="500">500</option>
            </Form.Select>
            <Button variant="primary" type="submit">
              Submit
            </Button>
          </Form>
        </Col>
        <Col>
          <Network 
            data={networkData} 
            width={networkDims.width}
            height={networkDims.height}
          />
          <ResultsTable data={articleMatches} />
        </Col>
      </Row>
    </Container>
  );
}

export default App;
