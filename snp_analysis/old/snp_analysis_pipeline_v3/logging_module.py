# logging_module.py
def log_message(message, level="normal", normal_part=""):
    color_map = {
        "info": ("\033[95m", "ℹ️ "),       # Magenta + Info Emoji
        "warning": ("\033[93m", "⚠️ "),    # Yellow + Warning Emoji
        "error": ("\033[91m", "❌ "),      # Red + Error Emoji
        "success": ("\033[92m", "✅ "),    # Green + Success Emoji
        "normal": ("\033[0m", ""),        # Default Reset (No Emoji)
        "blue": ("\033[94m", "🔵 ")       # Blue + Generic Emoji
    }
    
    color, emoji = color_map.get(level, ("\033[0m", ""))
    if normal_part:
        print(f"{color}{emoji}{message}\033[0m{normal_part}")
    else:
        print(f"{color}{emoji}{message}\033[0m")
