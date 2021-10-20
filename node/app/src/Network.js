import PropTypes from 'prop-types';
import ForceGraph2D from "react-force-graph-2d";

const Network = (props) => {
  return(
    <ForceGraph2D
      graphData={props.data}
      linkColor={() => '#eeeeee'}
      nodeLabel={'label'}
      nodeColor={'color'}
      nodeRelSize={8}
      width={props.width}
      height={props.height}
    />
  )
}

Network.propTypes = {
  data: PropTypes.object,
  height: PropTypes.number,
  width: PropTypes.number,
};

export default Network;
