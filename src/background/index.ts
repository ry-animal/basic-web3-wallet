// Service Worker background script
chrome.runtime.onInstalled.addListener(() => {
  console.log('Extension installed')
})

// Keep service worker active
chrome.runtime.onMessage.addListener((_message, _sender, _sendResponse) => {
  // Prefix unused parameters with underscore
  return true
})

export {} 