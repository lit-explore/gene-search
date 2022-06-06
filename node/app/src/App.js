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
  const [clusterInfo, setClusterInfo] = useState({})
  const [excludedGenes, setExcludedGenes] = useState("")
  const [networkData, setNetworkData] = useState({nodes: [], links: []})
  const [networkDims, setNetworkDims] = useState({width: 640, height: 480})

  const windowSize = useWindowSize();

  // set initial network dimensions, when windowSize hook is ready..
  useEffect(() => {
    if ((windowSize.height === undefined) || (windowSize.width === undefined)) {
      return
    }

    setNetworkDims({
      width: windowSize.width * 0.65,
      height: windowSize.height * 0.5
    })
  }, [windowSize.height, windowSize.width])

  function handleSubmit(event) {
    // construct query URL
    const geneInput = document.getElementById("input-genes").value 
    const genes = geneInput.trim().replace(/[,\n]/g, ",")

    const keyType = document.getElementById("gene-id-type").value 

    const pvalInput = document.getElementById("input-pvalues").value 
    const pvalues = pvalInput.trim().replace(/[,\n]/g, ",")
    
    const diseaseInput = document.getElementById("input-disease-mesh").value

    const maxArticles = document.getElementById("max-articles").value

    const baseURL = "http://lit.biodat.io/api/"
    let queryURL = `${baseURL}?genes=${genes}&keyType=${keyType}&max_articles=${maxArticles}`

    if (pvalues.length > 0) {
      queryURL = queryURL + `&pvalues=${pvalues}`
    }

    if (diseaseInput.length > 0) {
      queryURL = queryURL + `&disease=${diseaseInput}`
    }

    console.log(queryURL)

    // submit query to API
    axios.post(queryURL).then(resp => {
      console.log(resp)

      setArticleMatches(resp.data.network.nodes)
      setClusterInfo(resp.data.clusters)
      setExcludedGenes(resp.data.info.not_present_in_top_articles)
      setNetworkData(resp.data.network)

      // network tooltip text (title + genes)
      resp.data.network.nodes.forEach(entry => {
        entry.label = entry.citation + "<br /><b>" + entry.genes + "</b>"
      })
    }).catch(function (error) {
      console.log("[Error] Error encountered! See HTTP response in inspector for more info.");
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
            <Form.Label>Gene IDs</Form.Label>
            <Form.Control id="input-genes" as="textarea" placeholder="Gene IDs" style={{ height: '100px' }} />
            <Form.Text className="text-muted">
              Comma- or newline-separated gene IDs or symbols.
            </Form.Text>
            <Form.Select id="gene-id-type">
              <option value='symbol'>Gene Symbols</option>
              <option value='ensgene'>Ensembl Gene IDs</option>
            </Form.Select>
            <Form.Label>Gene P-values (Optional)</Form.Label>
            <Form.Control id="input-pvalues" as="textarea" placeholder="Gene P-values" style={{ height: '100px' }} />
            <Form.Text className="text-muted">
              Comma- or newline-separated P-values.
            </Form.Text>
            <Form.Label>Disease MeSH term (optional)</Form.Label>
            <Form.Control id="input-disease-mesh" placeholder='ex. "MESH:D009101"' />
            <Form.Text className="text-muted">
              Optional disease filter. If specified, only articles associated with the specified disease MeSH term will be considered.
            </Form.Text>
            <Form.Label>Article limit</Form.Label>
            <Form.Select id="max-articles" defaultValue="50" aria-label="Article Limit">
              <option value="25">25</option>
              <option value="50">50</option>
              <option value="100">100</option>
              <option value="150">150</option>
              <option value="200">200</option>
              <option value="250">250</option>
              <option value="500">500</option>
            </Form.Select>
            <Button variant="primary" type="submit">
              Submit
            </Button>
          </Form>
          <hr />
          <div id='cluster-results'>
          { Object.entries(clusterInfo).map(([clustId, cluster]) => (
            <div key={"cluster-" + clustId}>
              <span style={{color: cluster.color}}>
               Cluster { clustId } (n={cluster.num_articles} articles)
            </span><br />
            <ul>
              Common Genes:
              { Object.entries(cluster.genes).slice(0, 10).map(([gene, geneCount]) => (
                <li key={"cluster-" + clustId + "-gene-" + gene}>{gene} ({geneCount})</li>
              ))}
            </ul>
            </div>
          ))}
          </div>
          <hr />
          { excludedGenes !== "" &&
            <div id='infrequent-genes'>
              The following genes were not frequently found among the high-scoring 
              articles, and have been excluded from part of the analysis:
              <br />
              { excludedGenes }
            </div>
          }
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
