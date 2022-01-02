import PropTypes from 'prop-types';
import ForceGraph2D from "react-force-graph-2d";

const Network = (props) => {
  return(
      //nodeVal={node => Math.log(node.score + 1)}
    <ForceGraph2D
      graphData={props.data}
      linkColor={() => '#eeeeee'}
      linkLabel={'num_shared'}
      linkWidth={link => link.num_shared / 6}
      nodeLabel={'label'}
      nodeColor={'color'}
      nodeRelSize={8}
      nodeVal={node => Math.log(node.score + 1)/ 4}
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
