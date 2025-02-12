import { ethers } from 'ethers';

class WalletController {
  private wallet: ethers.HDNodeWallet | null = null;
  
  async createWallet(password: string) {
    const wallet = ethers.Wallet.createRandom() as ethers.HDNodeWallet;
    const encryptedWallet = await wallet.encrypt(password);
    await chrome.storage.local.set({ wallet: encryptedWallet });
    this.wallet = wallet;
    return wallet.address;
  }

  async unlockWallet(password: string) {
    const { wallet: encryptedWallet } = await chrome.storage.local.get('wallet');
    if (!encryptedWallet) throw new Error('No wallet found');
    
    this.wallet = await ethers.Wallet.fromEncryptedJson(encryptedWallet, password) as ethers.HDNodeWallet;
    return this.wallet?.address;
  }

  async signTransaction(tx: any) {
    if (!this.wallet) throw new Error('Wallet not unlocked');
    return await this.wallet.signTransaction(tx);
  }
}

const walletController = new WalletController();

chrome.runtime.onMessage.addListener((request, sender, sendResponse) => {
  (async () => {
    try {
      switch (request.type) {
        case 'CREATE_WALLET':
          const address = await walletController.createWallet(request.password);
          sendResponse({ success: true, address });
          break;
          
        case 'UNLOCK_WALLET':
          const unlockedAddress = await walletController.unlockWallet(request.password);
          sendResponse({ success: true, address: unlockedAddress });
          break;
          
        case 'SIGN_TRANSACTION':
          const signedTx = await walletController.signTransaction(request.transaction);
          sendResponse({ success: true, signedTx });
          break;
      }
    } catch (error) {
      sendResponse({ success: false, error: (error as Error).message });
    }
  })();
  return true;
}); 