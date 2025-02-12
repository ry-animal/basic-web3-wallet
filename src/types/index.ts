export interface TransactionRequest {
  to: string;
  from: string;
  value?: string;
  data?: string;
  gasLimit?: string;
  gasPrice?: string;
  nonce?: number;
}

export interface WalletMessage {
  type: 'CREATE_WALLET' | 'UNLOCK_WALLET' | 'SIGN_TRANSACTION' | 'REQUEST_ACCOUNTS';
  password?: string;
  transaction?: TransactionRequest;
}

export interface Web3Request {
  method: string;
  params: TransactionRequest[];
} 