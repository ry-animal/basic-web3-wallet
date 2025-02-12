/******/ (() => { // webpackBootstrap
/******/ 	"use strict";
/******/ 	// The require scope
/******/ 	var __webpack_require__ = {};
/******/ 	
/************************************************************************/
/******/ 	/* webpack/runtime/make namespace object */
/******/ 	(() => {
/******/ 		// define __esModule on exports
/******/ 		__webpack_require__.r = (exports) => {
/******/ 			if(typeof Symbol !== 'undefined' && Symbol.toStringTag) {
/******/ 				Object.defineProperty(exports, Symbol.toStringTag, { value: 'Module' });
/******/ 			}
/******/ 			Object.defineProperty(exports, '__esModule', { value: true });
/******/ 		};
/******/ 	})();
/******/ 	
/************************************************************************/
var __webpack_exports__ = {};
/*!******************************!*\
  !*** ./src/contentScript.ts ***!
  \******************************/
__webpack_require__.r(__webpack_exports__);
// Inject Web3 provider
window.addEventListener('load', () => {
    const provider = {
        isMetaMask: true,
        request: async ({ method, params }) => {
            switch (method) {
                case 'eth_requestAccounts':
                    return chrome.runtime.sendMessage({ type: 'REQUEST_ACCOUNTS' });
                case 'eth_sendTransaction':
                    return chrome.runtime.sendMessage({
                        type: 'SIGN_TRANSACTION',
                        transaction: params[0]
                    });
                // Add other RPC methods as needed
                default:
                    throw new Error(`Method ${method} not implemented`);
            }
        }
    };
    window.ethereum = provider;
});


/******/ })()
;
//# sourceMappingURL=contentScript.js.map