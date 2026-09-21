import { StrictMode } from 'react'
import { createRoot } from 'react-dom/client'
import './index.css'
import App from './App.jsx'

const standaloneRoot = document.getElementById('root')
if (standaloneRoot) {
  createRoot(standaloneRoot).render(
    <StrictMode>
      <App experimentId={standaloneRoot.dataset.experimentId || '1'} />
    </StrictMode>,
  )
}

document.querySelectorAll('.weaver-experiment-flow').forEach((element) => {
  createRoot(element).render(
    <StrictMode>
      <App experimentId={element.dataset.experimentId} />
    </StrictMode>,
  )
})
