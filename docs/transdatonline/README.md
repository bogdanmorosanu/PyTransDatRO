# TransDatOnline Documentation

This directory contains technical documentation and specifications for **TransDatOnline**, the web-facing coordinate transformation service and interface for Romania's geodetic infrastructure.

## Documents

1. [01_legacy_java_architecture.md](file:///d:/ProiecteRealizate/pyTransDatRO/repo/PyTransDatRO/docs/transdatonline/01_legacy_java_architecture.md)
   - Detailed architectural breakdown of the legacy Java / GWT implementation in [`transdatonline/`](file:///d:/ProiecteRealizate/pyTransDatRO/repo/PyTransDatRO/transdatonline).
   - Covers servlets (`JSONCooOpService`, `CoordinateOperationServiceImpl`), cartographic calculation classes, coordinate serialization models, MySQL database tables, and differences between legacy and modern geodetic models.

2. [02_legacy_service_api_specification.md](file:///d:/ProiecteRealizate/pyTransDatRO/repo/PyTransDatRO/docs/transdatonline/02_legacy_service_api_specification.md)
   - Definitive API specification for the public coordinate transformation endpoint `/transdatonline/cooOpService`.
   - Complete documentation of HTTP GET and POST request schemas, coordinate order/units (radians vs meters), JSON response formats, warning conditions (`Out of grid`, `No data on grid`, `Invalid coordinate data`), and numerical verification vectors.

3. [03_fastapi_compatibility_adapter_plan.md](file:///d:/ProiecteRealizate/pyTransDatRO/repo/PyTransDatRO/docs/transdatonline/03_fastapi_compatibility_adapter_plan.md)
   - Implementation architecture for integrating the backward compatibility layer into the new FastAPI service (`api/routes/legacy.py`).
   - Details route aliasing, error mapping, telemetry integration (`source=2`), and verification testing.

4. [04_legacy_web_ui_specification.md](file:///d:/ProiecteRealizate/pyTransDatRO/repo/PyTransDatRO/docs/transdatonline/04_legacy_web_ui_specification.md)
   - Comprehensive UI layout, visual design, and component hierarchy of the legacy GWT web interface (`MainForm`, `CRSPanel`, `OptionsPanel`, `TitlePanel`, `CommandsPanel`).
   - Documents CSS styling (`TransDatOnline.css`), typography, colors, modal dialogs (`WaitPanel`, `ReportPanel`), and user workflows.

5. [05_legacy_client_functionality_and_validation.md](file:///d:/ProiecteRealizate/pyTransDatRO/repo/PyTransDatRO/docs/transdatonline/05_legacy_client_functionality_and_validation.md)
   - Exhaustive analysis of client-side logic, line and delimiter tokenization, DMS / decimal degrees / radians conversions (`Angle.java`), and 2D/3D dimension validation.
   - Details coordinate order swapping (`NE` vs `EN`), GWT-RPC asynchronous execution lifecycle, output formatting, and the three-tier diagnostic reporting mechanism (`CooOpReport`).

6. [06_modern_web_app_transition_analysis.md](file:///d:/ProiecteRealizate/pyTransDatRO/repo/PyTransDatRO/docs/transdatonline/06_modern_web_app_transition_analysis.md)
   - Gap analysis and transition roadmap bridging the legacy GWT web app and the modern PyTransDatRO Web App (`webapp/`).
   - Details features preserved (DMS parser, NE coordinate order, diagnostic report), legacy bugs and quirks resolved, and REST API mapping.
