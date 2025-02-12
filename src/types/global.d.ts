import { Web3Request } from '../types';

interface Window {
  ethereum: {
    isMetaMask: boolean;
    request: (args: Web3Request) => Promise<string>;
  };
} 