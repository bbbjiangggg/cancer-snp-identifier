# logging_module.py
def log_message(message, level="normal", normal_part=""):
    color_map = {
        "info": ("\033[95m", "ℹ️"),       # Magenta, Info emoji
        "warning": ("\033[93m", "⚠️"),    # Yellow, Warning emoji
        "error": ("\033[91m", "❌"),      # Red, Error emoji
        "success": ("\033[92m", "✅"),    # Green, Success emoji
        "normal": ("\033[0m", ""),        # Reset to default
        "blue": ("\033[94m", "🔵")        # Blue, General note
    }

    color, emoji = color_map.get(level, ("\033[0m", ""))
    
    if normal_part:
        print(f"{color}{emoji} {message}\033[0m{normal_part}")
    else:
        print(f"{color}{emoji} {message}\033[0m")
