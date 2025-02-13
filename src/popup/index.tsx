import { createRoot } from 'react-dom/client'
import { Web3Provider } from '@/components/web3/Web3Provider'

const Popup = () => {
    return (
        <Web3Provider>
            <div className="w-[350px] h-[600px] bg-black flex flex-col items-center justify-center font-bold text-primary-500">
                <h1>Web3 Wallet</h1>
            </div>
        </Web3Provider>
    )
}

const container = document.getElementById('root')
if (container) {
    const root = createRoot(container)
    root.render(<Popup />)
} 