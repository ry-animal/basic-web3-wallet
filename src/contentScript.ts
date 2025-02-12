declare global {
  interface Window {
    ethereum: {
      isMetaMask: boolean;
      request: (args: { method: string; params: any[] }) => Promise<any>;
    };
  }
}

// Inject Web3 provider
window.addEventListener('load', () => {
  const provider = {
    isMetaMask: true,
    request: async ({ method, params }: { method: string; params: any[] }) => {
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

export {};  // Make this a module 