import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'
import path from 'path'
import { type UserConfig } from 'vite'

export default defineConfig({
  plugins: [react()],
  resolve: {
    alias: {
      '@': path.resolve(__dirname, './src')
    }
  },
  build: {
    rollupOptions: {
      input: {
        popup: 'src/popup/index.tsx',
        background: 'src/background/index.ts',
        content: 'src/content/index.ts'
      },
      output: {
        entryFileNames: '[name].js'
      }
    }
  }
} as UserConfig) 