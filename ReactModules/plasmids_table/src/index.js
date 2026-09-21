import React from 'react';
import ReactDOM from 'react-dom';
import './index.css';
import App from './App';
import PrimerApp from './PrimerApp';
import reportWebVitals from './reportWebVitals';

const plasmidsRoot = document.getElementById('root-plasmids');
if (plasmidsRoot) {
  ReactDOM.render(
    <React.StrictMode>
      <App />
    </React.StrictMode>,
    plasmidsRoot
  );
}

const primersRoot = document.getElementById('root-primers');
if (primersRoot) {
  ReactDOM.render(
    <React.StrictMode>
      <PrimerApp />
    </React.StrictMode>,
    primersRoot
  );
}

// If you want to start measuring performance in your app, pass a function
// to log results (for example: reportWebVitals(console.log))
// or send to an analytics endpoint. Learn more: https://bit.ly/CRA-vitals
reportWebVitals();
