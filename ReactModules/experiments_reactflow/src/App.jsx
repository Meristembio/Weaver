import { useState, useCallback, useEffect, useMemo, useRef } from 'react';
import {
  ReactFlow,
  Controls,
  Background,
  applyEdgeChanges,
  applyNodeChanges,
  Handle,
  Position,
  ReactFlowProvider,
  useReactFlow,
  MiniMap,
  Panel
} from '@xyflow/react';

import '@xyflow/react/dist/style.css';
import './App.css';

const CHILD_HEIGHT = 76;
const NODE_X_OFFSET = 360;
const SIBLING_GAP = 24;
const ROOT_GAP = 48;
const ROOT_COMPACT_GAP = 96;

const statusLabel = {
  V: 'Validated',
  UC: 'Under construction',
  NV: 'Not validated',
  RS: 'Reference',
};

const levelLabel = (item) => item.level === null ? 'No level' : `L${item.level}`;

const isPendingAssembly = (item) => item.status !== 'V' && item.status !== 'RS' && item.level !== null;

const sortPlasmids = (items) => {
  return [...items].sort((a, b) => {
    if ((a.level ?? 999) !== (b.level ?? 999)) return (a.level ?? 999) - (b.level ?? 999);
    return (a.weaver_id ?? 0) - (b.weaver_id ?? 0);
  });
};

const sortPlasmidsByBuildPriority = (items) => {
  return [...items].sort((a, b) => {
    const aPending = isPendingAssembly(a);
    const bPending = isPendingAssembly(b);
    if (aPending !== bPending) return aPending ? -1 : 1;
    if ((a.level ?? 999) !== (b.level ?? 999)) return (a.level ?? 999) - (b.level ?? 999);
    if (a.status !== b.status) return a.status.localeCompare(b.status);
    return (a.weaver_id ?? 0) - (b.weaver_id ?? 0);
  });
};

const getStatLists = (experiment) => {
  const plasmids = sortPlasmids(experiment.plasmids);
  return {
    total: sortPlasmidsByBuildPriority(experiment.plasmids),
    validated: sortPlasmids(plasmids.filter((item) => item.status === 'V')),
    pending: sortPlasmids(plasmids.filter(isPendingAssembly)),
    ready_to_build: sortPlasmids(plasmids.filter((item) => isPendingAssembly(item) && item.ready_to_build)),
    blocked: sortPlasmids(plasmids.filter((item) => isPendingAssembly(item) && !item.ready_to_build)),
  };
};

const statConfig = [
  { key: 'total', label: 'total' },
  { key: 'validated', label: 'validated' },
  { key: 'pending', label: 'pending' },
  { key: 'ready_to_build', label: 'ready' },
  { key: 'blocked', label: 'blocked' },
];

const groupByLevel = (items) => {
  return items.reduce((groups, item) => {
    const key = levelLabel(item);
    if (!groups.has(key)) groups.set(key, []);
    groups.get(key).push(item);
    return groups;
  }, new Map());
};

const getNodeProperties = (item) => {
  const baseStyle = {
    padding: '9px 12px',
    borderRadius: '6px',
    borderWidth: '1.5px',
    fontSize: '11px',
    fontWeight: 'bold',
    color: '#fff',
    width: '168px',
    textAlign: 'center',
    boxSizing: 'border-box',
  };

  let background = '#6c757d';
  let borderColor = '#adb5bd';
  if (item.status === 'V') {
    background = '#019256';
    borderColor = '#8fd1b6';
  } else if (item.status === 'UC') {
    background = '#0d6efd';
    borderColor = '#8bbcff';
  } else if (item.status === 'NV') {
    background = '#6f42c1';
    borderColor = '#c5a8ff';
  } else if (item.status === 'RS') {
    background = '#495057';
    borderColor = '#adb5bd';
  }

  let logicalType = 'default';
  if (!item.parent || item.parent.length === 0) logicalType = 'input';
  else if (!item.parts || item.parts.length === 0) logicalType = 'output';

  return {
    type: 'weaverNode',
    data: {
      logicalType,
      uuid: item.uuid,
      colony: item.colony,
      noColony: item.no_colony,
      status: item.status,
      level: item.level,
      name: item.name,
      weaverId: item.weaver_id,
      url: item.url,
      ligationRaw: item.ligation_raw,
      readyToBuild: item.ready_to_build,
    },
    style: { ...baseStyle, background, borderColor },
  };
};

const copyText = async (text) => {
  if (!text) return false;
  if (navigator.clipboard?.writeText && window.isSecureContext) {
    await navigator.clipboard.writeText(text);
    return true;
  }

  const textarea = document.createElement('textarea');
  textarea.value = text;
  textarea.setAttribute('readonly', '');
  textarea.style.position = 'fixed';
  textarea.style.top = '-1000px';
  textarea.style.left = '-1000px';
  document.body.appendChild(textarea);
  textarea.focus();
  textarea.select();

  let copied = false;
  try {
    copied = document.execCommand('copy');
  } finally {
    document.body.removeChild(textarea);
  }
  return copied;
};

const WeaverNode = ({ data, style }) => {
  const [copied, setCopied] = useState(false);
  const isInput = data.logicalType === 'input';
  const isOutput = data.logicalType === 'output';
  const isVerified = data.status === 'V';
  const colonyText = data.noColony ? 'NC' : (data.colony ?? 'NS');

  const handleBaseStyle = {
    width: '6px',
    height: '6px',
    background: '#555',
    border: '1px solid #fff',
    zIndex: 2,
  };

  const copyLigation = async (event) => {
    event.preventDefault();
    event.stopPropagation();
    const didCopy = await copyText(data.ligationRaw);
    if (!didCopy) return;
    setCopied(true);
    window.setTimeout(() => setCopied(false), 1200);
  };

  return (
    <div className="weaver-flow-node" style={style} title={statusLabel[data.status] || data.status}>
      {!isInput && (
        <Handle
          type="target"
          position={Position.Left}
          style={{ ...handleBaseStyle, left: '-12px' }}
          isConnectable={false}
        />
      )}

      {data.url && (
        <a
          href={data.url}
          target="_blank"
          rel="noopener noreferrer"
          className="weaver-flow-ref-button"
          onClick={(event) => event.stopPropagation()}
          title="Open plasmid"
        >
          <i className="bi bi-box-arrow-up-right" aria-hidden="true"></i>
        </a>
      )}

      {data.ligationRaw && data.status !== 'V' && data.level !== 0 && (
        <button
          type="button"
          className={`weaver-flow-copy-button${copied ? ' is-copied' : ''}`}
          onClick={copyLigation}
          onPointerDown={(event) => event.stopPropagation()}
          title={copied ? 'Copied' : 'Copy ligation data'}
        >
          <i className={`bi bi-${copied ? 'check2' : 'clipboard'}`} aria-hidden="true"></i>
        </button>
      )}

      <div className="weaver-flow-label">
        <div className="weaver-flow-node-meta">{levelLabel(data)} · {data.weaverId}</div>
        <div className="weaver-flow-node-name">{data.name}</div>
      </div>

      {isVerified && (
        <div className="weaver-flow-colony-badge" title={data.noColony ? 'No colony' : 'Working colony'}>
          {colonyText}
        </div>
      )}

      {!isOutput && (
        <Handle
          type="source"
          position={Position.Right}
          style={{ ...handleBaseStyle, right: '-12px' }}
          isConnectable={false}
        />
      )}
    </div>
  );
};

const animateNodePositions = (startNodes, targetNodes, setNodes, onComplete) => {
  const duration = 360;
  const startTime = performance.now();
  const startMap = new Map(startNodes.map((node) => [node.id, node]));
  const targetMap = new Map(targetNodes.map((node) => [node.id, node]));
  const allIds = new Set([...startMap.keys(), ...targetMap.keys()]);

  const frame = (currentTime) => {
    const elapsed = currentTime - startTime;
    const progress = Math.min(elapsed / duration, 1);
    const ease = 1 - Math.pow(1 - progress, 3);

    const currentFrameNodes = Array.from(allIds)
      .map((id) => {
        const startNode = startMap.get(id);
        const targetNode = targetMap.get(id);
        if (startNode && targetNode) {
          return {
            ...targetNode,
            position: {
              x: startNode.position.x + (targetNode.position.x - startNode.position.x) * ease,
              y: startNode.position.y + (targetNode.position.y - startNode.position.y) * ease,
            },
          };
        }
        return targetNode;
      })
      .filter(Boolean);

    setNodes(currentFrameNodes);
    if (progress < 1) requestAnimationFrame(frame);
    else if (onComplete) onComplete();
  };

  requestAnimationFrame(frame);
};

const nodeTypes = { weaverNode: WeaverNode };

const findRouteToPlasmid = (currentId, targetId, plasmidsByIdx, visited = new Set()) => {
  if (currentId === targetId) return [currentId];
  if (visited.has(currentId)) return null;
  visited.add(currentId);

  const current = plasmidsByIdx.get(currentId);
  if (!current?.parts?.length) return null;

  for (const childId of current.parts) {
    const childPath = findRouteToPlasmid(childId, targetId, plasmidsByIdx, new Set(visited));
    if (childPath) return [currentId, ...childPath];
  }

  return null;
};

const makeNodeFromItem = (item, id, position, isFocused = false) => {
  const props = getNodeProperties(item);
  const isReady = item.ready_to_build && item.status !== 'V';
  return {
    id,
    position,
    ...props,
    style: {
      ...props.style,
      borderColor: isFocused || isReady ? '#ffea00' : props.style.borderColor,
      borderWidth: isFocused || isReady ? '3px' : props.style.borderWidth,
      boxShadow: isReady ? '0 0 0 4px rgba(255, 234, 0, 0.72), 0 0 16px rgba(255, 234, 0, 0.42)' : undefined,
    },
    data: {
      ...props.data,
      label: `${item.weaver_id} | ${item.name}`,
      originalId: item.weaver_id,
    },
  };
};

const nodeIdForPath = (path) => path.join('-');

const expandedIdsForRoute = (route) => {
  const expanded = new Set();
  route.slice(0, -1).forEach((_, index) => {
    expanded.add(nodeIdForPath(route.slice(0, index + 1)));
  });
  return expanded;
};

const buildVisibleGraph = (experiment, plasmidsByIdx, expandedNodeIds, focusedPlasmidId, autoExpandBlocked) => {
  const route = focusedPlasmidId ? experiment.root_ids
    .map((rootId) => findRouteToPlasmid(rootId, focusedPlasmidId, plasmidsByIdx))
    .find(Boolean) : null;
  const focusedRootId = route?.[0] || null;
  const routeSet = new Set(route || []);
  const routeEdges = new Set();
  if (route) {
    route.slice(1).forEach((_, index) => {
      routeEdges.add(`${route[index]}-${route[index + 1]}`);
    });
  }

  const nodes = [];
  const edges = [];

  const shouldExpandNode = (item, nodeId, path) => {
    const belongsToFocusedRoot = !focusedRootId || path[0] === focusedRootId;
    return expandedNodeIds.has(nodeId) || (
      autoExpandBlocked && belongsToFocusedRoot && isPendingAssembly(item) && !item.ready_to_build
    );
  };

  const measureSubtree = (plasmidId, path, visited = new Set()) => {
    const item = plasmidsByIdx.get(plasmidId);
    const nodeId = nodeIdForPath(path);
    if (!item || visited.has(plasmidId)) return CHILD_HEIGHT;

    const childIds = shouldExpandNode(item, nodeId, path) ? item.parts || [] : [];
    if (!childIds.length) return CHILD_HEIGHT;

    const childrenHeight = childIds.reduce((height, childId, index) => {
      const childPath = [...path, childId];
      const childHeight = measureSubtree(childId, childPath, new Set([...visited, plasmidId]));
      return height + childHeight + (index < childIds.length - 1 ? SIBLING_GAP : 0);
    }, 0);

    return Math.max(CHILD_HEIGHT, childrenHeight);
  };

  const layoutSubtree = (plasmidId, path, depth, topY, visited = new Set(), forcedCenterY = null) => {
    const item = plasmidsByIdx.get(plasmidId);
    const nodeId = nodeIdForPath(path);
    if (!item || visited.has(plasmidId)) {
      return { height: CHILD_HEIGHT, centerY: topY + (CHILD_HEIGHT / 2) };
    }

    const shouldExpand = shouldExpandNode(item, nodeId, path);
    const childIds = shouldExpand ? item.parts || [] : [];
    const childLayouts = [];
    const measuredHeight = measureSubtree(plasmidId, path, visited);
    const subtreeTopY = forcedCenterY === null ? topY : forcedCenterY - (measuredHeight / 2);
    let cursorY = subtreeTopY;

    childIds.forEach((childId) => {
      const childPath = [...path, childId];
      const childLayout = layoutSubtree(
        childId,
        childPath,
        depth + 1,
        cursorY,
        new Set([...visited, plasmidId]),
      );
      childLayouts.push({ id: childId, path: childPath, layout: childLayout });
      cursorY += childLayout.height + SIBLING_GAP;
    });

    const contentHeight = childLayouts.length
      ? cursorY - subtreeTopY - SIBLING_GAP
      : CHILD_HEIGHT;
    const centerY = forcedCenterY !== null ? forcedCenterY : childLayouts.length
      ? (childLayouts[0].layout.centerY + childLayouts[childLayouts.length - 1].layout.centerY) / 2
      : subtreeTopY + (CHILD_HEIGHT / 2);

    nodes.push(makeNodeFromItem(
      item,
      nodeId,
      { x: depth * NODE_X_OFFSET, y: centerY },
      plasmidId === focusedPlasmidId,
    ));

    childLayouts.forEach((child) => {
      const childNodeId = nodeIdForPath(child.path);
      const routeEdgeKey = `${plasmidId}-${child.id}`;
      const isRouteEdge = routeEdges.has(routeEdgeKey) || (routeSet.has(plasmidId) && routeSet.has(child.id));
      edges.push({
        id: `e-${nodeId}-${childNodeId}`,
        source: nodeId,
        target: childNodeId,
        animated: isRouteEdge,
        style: {
          stroke: isRouteEdge ? '#ffea00' : '#7f8b96',
          strokeWidth: isRouteEdge ? 2 : 1,
        },
      });
    });

    return {
      height: Math.max(CHILD_HEIGHT, contentHeight),
      centerY,
    };
  };

  const expandedRootCount = experiment.root_ids.filter((rootId) => {
    const root = plasmidsByIdx.get(rootId);
    return root && shouldExpandNode(root, String(rootId), [rootId]);
  }).length;
  const useCompactRoots = expandedRootCount <= 1;
  let rootTopY = 0;
  experiment.root_ids.forEach((rootId, index) => {
    if (useCompactRoots) {
      const rootCenterY = (index * ROOT_COMPACT_GAP) + (CHILD_HEIGHT / 2);
      layoutSubtree(rootId, [rootId], 0, rootCenterY - (CHILD_HEIGHT / 2), new Set(), rootCenterY);
      return;
    }

    const rootLayout = layoutSubtree(rootId, [rootId], 0, rootTopY);
    rootTopY += rootLayout.height + ROOT_GAP;
  });

  return { nodes, edges };
};

const expandedIdsForFocusedPlasmid = (experiment, plasmidsByIdx, focusedPlasmidId) => {
  if (!focusedPlasmidId) return new Set();
  const route = experiment.root_ids
    .map((rootId) => findRouteToPlasmid(rootId, focusedPlasmidId, plasmidsByIdx))
    .find(Boolean);

  return route ? expandedIdsForRoute(route) : new Set();
};

function FlowInner({ experiment, focusedPlasmidId, onClearFocus }) {
  const { fitView } = useReactFlow();
  const plasmidsByIdx = useMemo(() => {
    return new Map(experiment.plasmids.map((item) => [item.weaver_id, item]));
  }, [experiment]);

  const initialGraph = useMemo(() => {
    return buildVisibleGraph(experiment, plasmidsByIdx, new Set(), focusedPlasmidId, false);
  }, [experiment, focusedPlasmidId, plasmidsByIdx]);

  const [expandedNodeIds, setExpandedNodeIds] = useState(new Set());
  const [autoExpandBlocked, setAutoExpandBlocked] = useState(false);
  const [nodes, setNodes] = useState(initialGraph.nodes);
  const [edges, setEdges] = useState(initialGraph.edges);
  const previousNodesRef = useRef(initialGraph.nodes);

  const visibleGraph = useMemo(() => {
    return buildVisibleGraph(experiment, plasmidsByIdx, expandedNodeIds, focusedPlasmidId, autoExpandBlocked);
  }, [autoExpandBlocked, expandedNodeIds, experiment, focusedPlasmidId, plasmidsByIdx]);

  useEffect(() => {
    const timer = window.setTimeout(() => {
      animateNodePositions(previousNodesRef.current, visibleGraph.nodes, setNodes, () => fitView({ duration: 650, padding: 0.22 }));
      previousNodesRef.current = visibleGraph.nodes;
      setEdges(visibleGraph.edges);
    }, 20);

    return () => window.clearTimeout(timer);
  }, [fitView, visibleGraph]);

  useEffect(() => {
    if (!focusedPlasmidId) return undefined;
    const timer = window.setTimeout(() => {
      setExpandedNodeIds(expandedIdsForFocusedPlasmid(experiment, plasmidsByIdx, focusedPlasmidId));
      setAutoExpandBlocked(true);
    }, 20);

    return () => window.clearTimeout(timer);
  }, [experiment, focusedPlasmidId, plasmidsByIdx]);

  useEffect(() => {
    const refit = () => {
      window.setTimeout(() => fitView({ duration: 300, padding: 0.18 }), 80);
    };
    document.addEventListener('shown.bs.collapse', refit);
    window.addEventListener('resize', refit);
    return () => {
      document.removeEventListener('shown.bs.collapse', refit);
      window.removeEventListener('resize', refit);
    };
  }, [fitView]);

  const onNodeClick = useCallback((event, clickedNode) => {
    const lookupId = clickedNode.data.originalId;
    const sourceItem = plasmidsByIdx.get(lookupId);
    if (!sourceItem?.parts?.length) return;

    setAutoExpandBlocked(false);
    setExpandedNodeIds((current) => {
      const next = new Set(current);
      if (next.has(clickedNode.id)) next.delete(clickedNode.id);
      else next.add(clickedNode.id);
      return next;
    });
  }, [plasmidsByIdx]);

  const collapseAll = useCallback(() => {
    setExpandedNodeIds(new Set());
    setAutoExpandBlocked(false);
    onClearFocus();
    window.setTimeout(() => fitView({ duration: 650, padding: 0.22 }), 80);
  }, [fitView, onClearFocus]);

  return (
    <ReactFlow
      nodes={nodes}
      nodeTypes={nodeTypes}
      edges={edges}
      onNodesChange={(changes) => setNodes((currentNodes) => applyNodeChanges(changes, currentNodes))}
      onEdgesChange={(changes) => setEdges((currentEdges) => applyEdgeChanges(changes, currentEdges))}
      onNodeClick={onNodeClick}
      nodesDraggable={false}
      fitView
      colorMode="system"
    >
      <Panel position="top-left">
        <button
          type="button"
          className="weaver-flow-map-button"
          onClick={collapseAll}
          title="Collapse all"
        >
          <i className="bi bi-arrows-collapse" aria-hidden="true"></i>
          <span>Collapse</span>
        </button>
      </Panel>
      <Background />
      <Controls />
      <MiniMap
        position="bottom-right"
        zoomable
        pannable
        nodeStrokeColor={(node) => node.style?.background || '#888'}
        nodeColor={(node) => node.style?.background || '#eee'}
        nodeBorderRadius={2}
      />
    </ReactFlow>
  );
}

function ExperimentFlow({ experimentId }) {
  const [state, setState] = useState({ loading: true, experiment: null, error: '' });
  const [activeStat, setActiveStat] = useState('');
  const [focusedPlasmidId, setFocusedPlasmidId] = useState(null);

  useEffect(() => {
    let cancelled = false;
    fetch('/inventory/api/experiments-map/')
      .then((response) => {
        if (!response.ok) throw new Error(`HTTP ${response.status}`);
        return response.json();
      })
      .then((payload) => {
        const experiments = payload.projects.flatMap((project) => project.experiments);
        const experiment = experiments.find((item) => String(item.id) === String(experimentId));
        if (!experiment) throw new Error('Experiment data not found');
        if (!cancelled) setState({ loading: false, experiment, error: '' });
      })
      .catch((error) => {
        if (!cancelled) setState({ loading: false, experiment: null, error: error.message });
      });

    return () => {
      cancelled = true;
    };
  }, [experimentId]);

  if (state.loading) return <div className="weaver-flow-state">Loading map...</div>;
  if (state.error) return <div className="weaver-flow-state">Unable to load map: {state.error}</div>;
  if (!state.experiment?.plasmids?.length) return <div className="weaver-flow-state">No plasmids to map.</div>;

  const { stats } = state.experiment;
  const statLists = getStatLists(state.experiment);
  const selectedStat = statConfig.find((item) => item.key === activeStat);
  const selectedPlasmids = activeStat ? statLists[activeStat] : [];
  const selectedGroups = groupByLevel(selectedPlasmids);

  return (
    <div className={`weaver-flow-shell${state.experiment.archived ? ' is-archived' : ''}`}>
      <div className="weaver-flow-stats" aria-label="Experiment assembly stats">
        {statConfig.map((stat) => (
          <button
            type="button"
            key={stat.key}
            className={`weaver-flow-stat-chip${activeStat === stat.key ? ' is-active' : ''}`}
            onClick={() => setActiveStat((current) => current === stat.key ? '' : stat.key)}
            aria-expanded={activeStat === stat.key}
          >
            <strong>{stats[stat.key]}</strong> {stat.label}
          </button>
        ))}
        <span><strong>{stats.progress}%</strong> complete</span>
        {state.experiment.archived && (
          <span className="weaver-flow-archived-badge">
            <i className="bi bi-archive" aria-hidden="true"></i> Archived
          </span>
        )}
      </div>
      {selectedStat && (
        <div className="weaver-flow-stat-panel">
          <div className="weaver-flow-stat-panel-header">
            <strong>{selectedPlasmids.length}</strong>
            <span>{selectedStat.label}</span>
            <button
              type="button"
              className="weaver-flow-stat-panel-close"
              onClick={() => setActiveStat('')}
              title="Close list"
            >
              <i className="bi bi-x-lg" aria-hidden="true"></i>
            </button>
          </div>
          <div className="weaver-flow-stat-list">
            {selectedPlasmids.length ? Array.from(selectedGroups.entries()).map(([level, plasmids]) => (
              <div className="weaver-flow-stat-group" key={`${selectedStat.key}-${level}`}>
                <div className="weaver-flow-stat-group-title">{level}</div>
                {plasmids.map((plasmid) => (
                  <button
                    type="button"
                    key={`${selectedStat.key}-${plasmid.uuid}`}
                    className={`weaver-flow-stat-row${focusedPlasmidId === plasmid.weaver_id ? ' is-selected' : ''}`}
                    onClick={() => {
                      setFocusedPlasmidId(plasmid.weaver_id);
                      setActiveStat('');
                    }}
                    title="Show route in map"
                  >
                    <span className={`weaver-flow-stat-status status-${plasmid.status.toLowerCase()}`}>
                      {plasmid.status}
                    </span>
                    <span className="weaver-flow-stat-name">{plasmid.name}</span>
                    <strong className="weaver-flow-stat-id">{plasmid.weaver_id}</strong>
                    {(plasmid.colony !== null || plasmid.no_colony) && (
                      <span className="weaver-flow-stat-colony">
                        {plasmid.no_colony ? 'NC' : `c ${plasmid.colony}`}
                      </span>
                    )}
                  </button>
                ))}
              </div>
            )) : (
              <div className="weaver-flow-stat-empty">No plasmids in this group.</div>
            )}
          </div>
        </div>
      )}
      <div className="weaver-flow-canvas">
        <ReactFlowProvider>
          <FlowInner
            experiment={state.experiment}
            focusedPlasmidId={focusedPlasmidId}
            onClearFocus={() => setFocusedPlasmidId(null)}
          />
        </ReactFlowProvider>
      </div>
    </div>
  );
}

export default ExperimentFlow;
