class PopupUI {
  async initialize() {
    console.log('PopupUI initializing...');  // Debug log
    
    const createButton = document.getElementById('createWallet');
    const unlockButton = document.getElementById('unlockWallet');
    
    if (!createButton || !unlockButton) {
      console.error('Could not find buttons!');  // Debug log
      return;
    }

    createButton.addEventListener('click', () => {
      console.log('Create wallet clicked');  // Debug log
      const password = prompt('Enter a password for your wallet:');
      if (password) {
        this.createWallet(password);
      }
    });

    unlockButton.addEventListener('click', () => {
      console.log('Unlock wallet clicked');  // Debug log
      const password = prompt('Enter your wallet password:');
      if (password) {
        this.unlockWallet(password);
      }
    });
  }

  private async createWallet(password: string) {
    try {
      const response = await chrome.runtime.sendMessage({ 
        type: 'CREATE_WALLET', 
        password 
      });
      
      if (response.success) {
        this.showWalletView(response.address);
      } else {
        alert('Failed to create wallet: ' + response.error);
      }
    } catch (error) {
      console.error('Create wallet error:', error);  // Debug log
      alert('Error creating wallet');
    }
  }

  private async unlockWallet(password: string) {
    try {
      const response = await chrome.runtime.sendMessage({ 
        type: 'UNLOCK_WALLET', 
        password 
      });
      
      if (response.success) {
        this.showWalletView(response.address);
      } else {
        alert('Failed to unlock wallet: ' + response.error);
      }
    } catch (error) {
      console.error('Unlock wallet error:', error);  // Debug log
      alert('Error unlocking wallet');
    }
  }

  private showWalletView(address: string) {
    console.log('Showing wallet view for address:', address);  // Debug log
    const welcomeScreen = document.getElementById('welcome');
    const walletScreen = document.getElementById('walletView');
    const addressElement = document.getElementById('address');

    if (!welcomeScreen || !walletScreen || !addressElement) {
      console.error('Missing UI elements');  // Debug log
      return;
    }

    welcomeScreen.classList.add('hidden');
    walletScreen.classList.remove('hidden');
    addressElement.textContent = address;
  }
}

// Wait for DOM to be ready
document.addEventListener('DOMContentLoaded', () => {
  console.log('DOM loaded, initializing PopupUI');  // Debug log
  new PopupUI().initialize();
}); 