# Experiments ReactFlow

React/Vite module for the experiments assembly map shown in Django at
`/inventory/experiments`.

## Build

From this directory:

```sh
npm install
npm run build
```

The Vite build writes the compiled assets directly to:

```text
../../Django/static/experiments-reactflow
```

Django loads those generated files from
`inventory/templates/inventory/parts/base.html` when the request path is
`/inventory/experiments`.
