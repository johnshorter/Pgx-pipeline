"""Shared Jinja2 + WeasyPrint helpers for both report tiers."""

import logging
from pathlib import Path

from jinja2 import Environment, FileSystemLoader

from config.settings import TEMPLATES_DIR

logger = logging.getLogger(__name__)


def make_environment() -> Environment:
    return Environment(
        loader=FileSystemLoader(str(TEMPLATES_DIR)),
        autoescape=True,
        trim_blocks=True,
        lstrip_blocks=True,
    )


def render_to_html(template_name: str, context: dict, output_path: Path) -> Path:
    env = make_environment()
    html = env.get_template(template_name).render(**context)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html, encoding="utf-8")
    return output_path


def render_pdf_if_possible(html: str, output_path: Path) -> Path | None:
    """Render PDF via WeasyPrint. Returns the output path on success, None
    if WeasyPrint is unavailable, and writes a placeholder .txt next to the
    expected PDF in that case."""
    try:
        from weasyprint import HTML
        HTML(string=html).write_pdf(str(output_path))
        return output_path
    except ImportError:
        _write_placeholder(output_path)
        return None
    except Exception as e:
        logger.warning("PDF generation failed: %s", e)
        _write_placeholder(output_path)
        return None


def _write_placeholder(pdf_path: Path) -> None:
    placeholder = pdf_path.with_name(pdf_path.stem + ".pdf_not_generated.txt")
    placeholder.write_text(
        "PDF generation requires WeasyPrint.\n"
        "Install it with: pip install weasyprint\n"
        "Then re-run the report generation.\n",
        encoding="utf-8",
    )


def html_and_pdf(
    template_name: str,
    context: dict,
    output_dir: Path,
    name: str,
) -> dict[str, Path]:
    """Render template to HTML and (best-effort) PDF. Returns paths."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    env = make_environment()
    html = env.get_template(template_name).render(**context)

    html_path = output_dir / f"{name}.html"
    html_path.write_text(html, encoding="utf-8")

    pdf_path = output_dir / f"{name}.pdf"
    pdf_result = render_pdf_if_possible(html, pdf_path)

    out: dict[str, Path] = {"html": html_path}
    if pdf_result is not None:
        out["pdf"] = pdf_result
    return out
