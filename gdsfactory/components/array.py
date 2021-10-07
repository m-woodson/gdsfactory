from typing import Tuple

import gdsfactory as gf
from gdsfactory.cell import cell
from gdsfactory.component import Component
from gdsfactory.components.straight import straight


@cell
def array(
    component: gf.types.ComponentOrFactory = straight,
    spacing: Tuple[float, float] = (150.0, 150.0),
    columns: int = 6,
    rows: int = 1,
) -> Component:
    """Returns an array of components.

    Args:
        component: to replicate
        spacing: x, y spacing
        columns: in x
        rows: in y


    .. code::

        2 rows x 4 columns
         ___        ___       ___          ___
        |   |      |   |     |   |        |   |
        |___|      |___|     |___|        |___|


         ___        ___       ___          ___
        |   |      |   |     |   |        |   |
        |___|      |___|     |___|        |___|


    """
    if rows > 1 and spacing[1] == 0:
        raise ValueError(f"rows = {rows} > 1 require spacing[1] > 0")

    if columns > 1 and spacing[0] == 0:
        raise ValueError(f"columns = {columns} > 1 require spacing[0] > 0")

    c = Component()
    component = component() if callable(component) else component
    c.add_array(component, columns=columns, rows=rows, spacing=spacing)

    for col in range(columns):
        for row in range(rows):
            for port in component.ports.values():
                name = f"{port.name}_{row+1}_{col+1}"
                c.add_port(name, port=port)
                c.ports[name].move((col * spacing[0], row * spacing[1]))
    return c


if __name__ == "__main__":

    c2 = array(rows=2, columns=2, spacing=(100, 100))
    c2.show(show_ports=True)
