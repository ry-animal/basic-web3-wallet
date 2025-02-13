import type { Config } from 'tailwindcss'

export default {
  content: [
    "./src/**/*.{js,jsx,ts,tsx}",
    "./public/popup.html"
  ],
  theme: {
    extend: {
      colors: {
        primary: {
          DEFAULT: '#3B82F6',
          dark: '#2563EB'
        }
      }
    }
  },
  plugins: []
} satisfies Config 