// Content script that runs on web pages
console.log('Content script loaded')

// Listen for messages from the extension
chrome.runtime.onMessage.addListener((_message, _sender, _sendResponse) => {
  // Prefix unused parameters with underscore
  return true
})

export {} 