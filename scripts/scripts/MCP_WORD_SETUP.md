# MCP Word Document Editing Setup Guide

## Quick Setup for Claude Desktop

### Method 1: Direct NPX Configuration (Simplest)

1. Open Claude Desktop settings
2. Navigate to MCP Server configuration
3. Add the following configuration:

```json
{
  "mcpServers": {
    "office-word": {
      "command": "npx",
      "args": ["-y", "@gongrzhe/office-word-mcp-server"],
      "env": {}
    }
  }
}
```

4. Restart Claude Desktop

### Method 2: Using Smithery (Requires API Key)

1. Get a free API key from https://smithery.ai/account/api-keys
2. Install via command line:
```bash
npx @smithery/cli install @GongRzhe/Office-Word-MCP-Server --client claude
```
3. Follow the prompts and enter your API key when asked

### Method 3: Local Installation (Advanced)

If you prefer a local installation, use the provided setup script:
```bash
./scripts/setup_word_mcp.sh
```

## Features Available

Once configured, you can use natural language commands to:

- **Create Documents**: "Create a new Word document called 'report.docx'"
- **Add Content**: "Add a heading 'Quarterly Report' and three paragraphs"
- **Format Text**: "Make the first paragraph bold and blue"
- **Insert Tables**: "Add a 4x4 table with sales data"
- **Add Images**: "Insert the company logo image"
- **Search & Replace**: "Replace all instances of 'draft' with 'final'"
- **Extract Content**: "Extract all text from the document"
- **Document Properties**: "Show document statistics and properties"

## Testing the Setup

After configuration, test with these commands in Claude:
1. "Create a test Word document"
2. "Add a title 'Test Document' to the document"
3. "Save the document as 'test.docx'"

## Troubleshooting

### If MCP server doesn't work:
1. Ensure Node.js is installed (version 16+)
2. Check Claude Desktop logs for errors
3. Verify the configuration is correctly added to Claude Desktop settings

### Common Issues:
- **"Command not found"**: Install Node.js first
- **"Permission denied"**: Run with appropriate permissions
- **"Module not found"**: Clear npm cache: `npm cache clean --force`

## Additional Resources

- Official MCP Server: https://github.com/GongRzhe/Office-Word-MCP-Server
- Smithery Package: https://server.smithery.ai/@GongRzhe/Office-Word-MCP-Server/mcp
- MCP Documentation: https://modelcontextprotocol.io/

## Summary

The Office-Word-MCP-Server enables Claude to directly manipulate Word documents through natural language commands. Choose the setup method that works best for your environment - the NPX configuration (Method 1) is the simplest and requires no installation.