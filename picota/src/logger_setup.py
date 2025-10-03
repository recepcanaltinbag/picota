import logging

def setup_logger_from_config(cfg):
    """
    Config üzerinden logger oluşturur.
    cfg.logging.log_file ve cfg.logging.level beklenir.
    """
    logger = logging.getLogger(cfg.logging.logger_name)
    if logger.handlers:
        return logger  # zaten oluşturulmuşsa tekrar ekleme

    # Level çevirisi string -> logging level
    level_str = getattr(cfg.logging, "level", "DEBUG")
    level = getattr(logging, level_str.upper(), logging.DEBUG)

    # Console handler
    ch = logging.StreamHandler()
    ch.setLevel(level)
    ch_formatter = logging.Formatter("[%(levelname)s] %(message)s")
    ch.setFormatter(ch_formatter)
    logger.addHandler(ch)

    # File handler opsiyonel
    log_file = getattr(cfg.logging, "log_file", None)
    if log_file:
        fh = logging.FileHandler(log_file)
        fh.setLevel(level)
        fh_formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
        fh.setFormatter(fh_formatter)
        logger.addHandler(fh)

    logger.setLevel(level)
    return logger

