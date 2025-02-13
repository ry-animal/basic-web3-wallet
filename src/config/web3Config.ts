import { createConfig, http } from 'wagmi'
import { mainnet, sepolia } from 'wagmi/chains'

export const config = createConfig({
  chains: [mainnet, sepolia],
  transports: {
    [mainnet.id]: http(),
    [sepolia.id]: http()
  }
})

export const SUPPORTED_CHAINS = [
  {
    id: mainnet.id,
    name: 'Ethereum',
    rpcUrl: 'https://eth-mainnet.g.alchemy.com/v2/YOUR_API_KEY'
  },
  {
    id: sepolia.id,
    name: 'Sepolia',
    rpcUrl: 'https://eth-sepolia.g.alchemy.com/v2/YOUR_API_KEY'
  }
] 