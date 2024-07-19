"""Based on Gerber file spec.

https://www.ucamco.com/files/downloads/file_en/456/gerber-layer-format-specification-revision-2022-02_en.pdf.

See Also:
- https://github.com/opiopan/pcb-tools-extension
- https://github.com/jamesbowman/cuflow/blob/master/gerber.py
"""

from pathlib import Path
from typing import Literal

from pydantic import BaseModel, Field

from gdsfactory import Component
from gdsfactory.typings import ComponentSpec, LayerSpec


class GerberLayer(BaseModel):
    name: str
    function: list[str]
    polarity: Literal["Positive", "Negative"]


# todo we should write this initially for mm only and then add functionality for inches
# todo we need to address using different resolutions and int sizes
class GerberOptions(BaseModel):
    header: list[str] | None = None
    mode: Literal["mm", "in"] = "mm"
    resolution: float = 1e-6
    int_size: int = 4


# For generating a gerber job json file
class BoardOptions(BaseModel):
    size: tuple | None = (None,)
    n_layers: int = (2,)


resolutions = {1e-3: 3, 1e-4: 4, 1e-5: 5, 1e-6: 6}


def number(n) -> str:
    i = int(round(n * 10000))
    return "%07d" % i


# D02 moves the "cursor" to the point, "D01" draws it
def points(pp: list):
    out = ""
    d = "D02"
    for x, y in pp:
        out += f"X{number(x)}Y{number(y)}{d}*\n"
        d = "D01"
    return out


# we are reserving D10 (custom aperture) for drawing polygons
def rect(x0, y0, x1, y1):
    return "D10*\n" + points([(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)])


def linestring(pp):
    return "D10*\n" + points(pp)


def polygon(pp):
    return "G36*\n" + points(pp) + "G37*\n" + "\n"


def circle(center_x, center_y):
    # requires that the aperture has already been selected
    return f"X{center_x}Y{center_y}D03*\n"


def get_circle_components(
    component: ComponentSpec, circles: dict, origin_offset: tuple[float, float]
) -> dict:
    """Recursively searches for and collects information about "circle" components within a given component.

    Args:
        component: The component to search within.
        circles: A dictionary to store the found circle components and the locations of their shifted centers.
        origin_offset: A tuple representing the offset of the current component's origin relative to the
        top-level origin.

    Returns:
        A dictionary containing the found circle components (as keys) and their shifted centers (as values).

    This function traverses the hierarchy of components, looking for components whose names start with "circle."
    For each such component, it calculates the shifted center (taking into account the origin offsets) and stores this
    information in the 'circles' dict.

    If a component contains other components, the function recursively calls itself on those nested components to
    continue the search.
    """
    for ref in component.insts:
        ref_offset = (
            origin_offset[0] + ref.dcenter.x,
            origin_offset[1] + ref.dcenter.y,
        )

        if ref.cell.name.startswith("circle"):
            center = ref.dcenter
            shifted_center = (center.x + origin_offset[0], center.y + origin_offset[1])

            print(f"{ref.cell.name}, relative center=({center})")

            if ref.cell not in circles.keys():
                circles[ref.cell] = [shifted_center]
            else:
                circles[ref.cell].append(shifted_center)

        elif len(ref.cell.insts) > 0:
            get_circle_components(ref.cell, circles, origin_offset=ref_offset)
    return circles


class Aperture(BaseModel):
    """Represents an aperture with a specified radius, layer, and locations.

    Attributes:
        radius: A float representing the radius of the aperture.
        layer: The layer on which the aperture is located. #todo
        locations: A list of tuples, where each tuple represents the (x, y) coordinates of an aperture center.

    This class serves as a data model to store information about an aperture, including its size (radius), the layer it
    belongs to, and the specific locations (centers) where it exists.
    """

    radius: float
    layer: int  # LayerSpec todo
    locations: list[tuple]


def get_aperture_list(circles: dict[ComponentSpec, list[tuple]]) -> list:
    """Creates a list of Aperture objects from a dictionary of circle components.

    Args:
        circles: A dictionary where keys are ComponentSpec objects representing circle components, and values are lists
        of tuples representing the shifted centers of those circles.

    Returns:
        A list of Aperture objects created using the information from the 'circles' dictionary.

    This function iterates through the 'circles' dictionary. For each circle component, it extracts the radius from the
    component's settings, the layer information, and the list of shifted centers. It then creates an Aperture object
    using these values and adds it to the 'aperture_list'.

    The function assumes that the ComponentSpec objects have a 'settings' attribute with a "radius" key and a 'layer()'
    method to retrieve layer information.
    """
    aperture_list = []
    for k, v in circles.items():
        # todo do we need a way to handle multiple layers from a given component?
        ap = Aperture(radius=k.settings["radius"], layer=k.layer(), locations=v)
        aperture_list.append(ap)
    return aperture_list


def get_apertures_by_layer(
    aperture_list: list[Aperture],
) -> dict[LayerSpec, list[Aperture]]:
    """Organizes a list of Aperture objects into a dictionary based on their layers.

    Args:
        aperture_list: A list of Aperture objects.

    Returns:
        A dictionary where keys are LayerSpec objects representing the layers, and values are lists of Aperture objects
        belonging to those respective layers.

    For each Aperture in the 'aperture_list', this function extracts the layer information and adds the Aperture to the
    corresponding list in the 'layer_to_apertures' dictionary. If a layer is encountered for the first time, a new entry
    is created in the dictionary with that layer as the key and an initial list containing the current Aperture.

    The function effectively groups Apertures based on their layers.
    """
    layer_to_apertures = {}

    for aperture in aperture_list:
        layer = aperture.layer
        if layer in layer_to_apertures.keys():
            layer_to_apertures[layer].append(aperture)
        else:
            layer_to_apertures[layer] = [aperture]

    return layer_to_apertures


def to_gerber(
    top_component: Component,
    dirpath: Path,
    layermap_to_gerber_layer: dict[tuple[int, int], GerberLayer],
    options: GerberOptions = Field(default_factory=dict),
) -> None:
    """Writes each layer to a different Gerber file.

    Args:
        top_component: to export.
        dirpath: directory path.
        layermap_to_gerber_layer: map of GDS layer to GerberLayer.
        options: to save.
            header: List[str] | None = None
            mode: Literal["mm", "in"] = "mm"
            resolution: float = 1e-6
            int_size: int = 4
    """
    # # Each layer and a list of the polygons (as lists of points) on that layer
    # layer_to_polygons = top_component.get_polygons_points()

    # Gets a dictionary of all circles components (keys) and a list of tuples with their reference center locations (values)
    circles_dict = get_circle_components(
        component=top_component, circles={}, origin_offset=(0, 0)
    )

    for k, v in circles_dict.items():
        print(f"{k}: {v}")

    aperture_list = get_aperture_list(circles_dict)
    print(aperture_list)

    layer_to_apertures = get_apertures_by_layer(aperture_list)
    print(layer_to_apertures)

    # for layer_tup, layer in layermap_to_gerber_layer.items():
    #     filename = (dirpath / layer.name.replace(" ", "_")).with_suffix(".gbr")
    #
    #     with open(filename, "w+") as f:
    #         header = options.header or [
    #             "Gerber file generated by gdsfactory",
    #             f"Component: {component.name}",
    #         ]
    #
    #         # Write file spec info
    #         f.write("%TF.FileFunction," + ",".join(layer.function) + "*%\n")
    #         f.write(f"%TF.FilePolarity,{layer.polarity}*%\n")
    #
    #         digits = resolutions[options.resolution]
    #         f.write(f"%FSLA{options.int_size}{digits}Y{options.int_size}{digits}X*%\n")
    #
    #         # Write header comments
    #         f.writelines([f"G04 {line}*\n" for line in header])
    #
    #         # Setup units/mode
    #         units = options.mode.upper()
    #         f.write(f"%MO{units}*%\n")
    #         f.write("%LPD*%")
    #
    #         f.write("G01*\n")
    #
    #         # Aperture definition
    #         f.write("%ADD10C,0.050000*%\n")
    #
    #         # Only supports polygons and circles for now
    #         if layer_tup in layer_to_polygons.keys():
    #             for poly in layer_to_polygons[layer_tup.layer]:
    #                 f.write(polygon(poly))
    #
    #         if layer_tup in layer_to_circles.keys():
    #             for circle in layer_to_circles[layer_tup.layer]:
    #                 f.write(polygon(poly))
    #
    #         # File end
    #         f.write("M02*\n")


if __name__ == "__main__":
    import gdsfactory as gf
    from gdsfactory.config import PATH
    from gdsfactory.technology import (
        LayerMap,
        LayerView,
        LayerViews,
    )
    from gdsfactory.typings import Layer

    class LayerMapPCB(LayerMap):
        F_Cu: Layer = (1, 0)
        In1_Cu: Layer = (2, 0)
        In2_Cu: Layer = (3, 0)
        B_Cu: Layer = (4, 0)
        F_Silkscreen: Layer = (11, 0)
        F_Mask: Layer = (21, 0)
        B_Mask: Layer = (22, 0)
        Edge_Cuts: Layer = (31, 0)

        DEVREC: Layer = (68, 0)
        PORT: Layer = (1, 10)
        PORTE: Layer = (1, 11)

    LAYER = LayerMapPCB

    class PCBViews(LayerViews):
        F_Cu: LayerView = LayerView(
            name="F_Cu",
            layer=tuple(LAYER.F_Cu),
            color="red",
        )
        In1_Cu: LayerView = LayerView(
            name="In1_Cu",
            layer=tuple(LAYER.In1_Cu),
            color="limegreen",
        )
        In2_Cu: LayerView = LayerView(
            name="In2_Cu",
            layer=tuple(LAYER.In2_Cu),
            color="goldenrod",
        )
        B_Cu: LayerView = LayerView(
            name="B_Cu",
            layer=tuple(LAYER.B_Cu),
            color="blue",
        )
        F_Silkscreen: LayerView = LayerView(
            name="F_Silkscreen",
            layer=tuple(LAYER.F_Silkscreen),
            color="khaki",
        )
        F_Mask: LayerView = LayerView(
            name="F_Mask",
            layer=tuple(LAYER.F_Mask),
            color="violet",
        )
        B_Mask: LayerView = LayerView(
            name="B_Mask",
            layer=LAYER.B_Mask,
            color="aqua",
        )
        Edge_Cuts: LayerView = LayerView(
            name="Edge_Cuts",
            layer=LAYER.Edge_Cuts,
            color="gold",
        )

    LAYER_VIEWS = PCBViews()

    # def get_pcb_layer_stack(
    #     copper_thickness: float = 0.035,
    #     core_thickness: float = 1,
    # ):
    #     return LayerStack(
    #         layers=dict(
    #             top_cu=LayerLevel(
    #                 layer=LAYER.F_Cu,
    #                 thickness=copper_thickness,
    #                 zmin=0.0,
    #                 material="cu",
    #             ),
    #             inner1_cu=LayerLevel(
    #                 layer=LAYER.In1_Cu,
    #                 thickness=copper_thickness,
    #                 zmin=0.0,
    #                 material="cu",
    #             ),
    #             inner_core=LayerLevel(
    #                 layer=LAYER.Edge_Cuts,
    #                 thickness=core_thickness,
    #                 zmin=-core_thickness / 2,
    #                 material="fr4",
    #             ),
    #             inner2_cu=LayerLevel(
    #                 layer=LAYER.In2_Cu,
    #                 thickness=copper_thickness,
    #                 zmin=0.0,
    #                 material="cu",
    #             ),
    #             bottom_cu=LayerLevel(
    #                 layer=LAYER.B_Cu,
    #                 thickness=copper_thickness,
    #                 zmin=0.0,
    #                 material="cu",
    #             ),
    #         )
    #     )
    #
    # LAYER_STACK = get_pcb_layer_stack()

    layermap_to_gerber = {
        LAYER.F_Cu: GerberLayer(
            name="F_Cu", function=["Copper", "L1", "Top"], polarity="Positive"
        ),
        LAYER.B_Cu: GerberLayer(
            name="B_Cu", function=["Copper", "L2", "Bot"], polarity="Positive"
        ),
        LAYER.F_Silkscreen: GerberLayer(
            name="F_Silkscreen", function=["Legend", "Top"], polarity="Positive"
        ),
        LAYER.F_Mask: GerberLayer(
            name="F_Mask", function=["SolderMask", "Top"], polarity="Negative"
        ),
        LAYER.B_Mask: GerberLayer(
            name="B_Mask", function=["SolderMask", "Bot"], polarity="Negative"
        ),
        LAYER.Edge_Cuts: GerberLayer(
            name="Edge_Cuts", function=["Profile"], polarity="Positive"
        ),
    }

    #     from gdsfactory.install import install_klayout_technology
    #     from gdsfactory.technology.klayout_tech import KLayoutTechnology
    #     tech_dir = (Path(__file__) / "..").resolve() / "klayout"
    #     pcb_tech = KLayoutTechnology(name='PCB', layer_views=PCBViews(), layer_map=LAYER)
    #     #pcb_tech.technology.dbu = 1e-3
    #     pcb_tech.write_tech(tech_dir=str(tech_dir))

    # install_klayout_technology(tech_dir=tech_dir, tech_name="PCB")

    c = gf.Component(name="top")
    c << gf.components.circle(layer=LAYER.F_Cu, radius=10)
    c << gf.components.rectangle(layer=LAYER.F_Cu)

    d = gf.Component(name="dummy_component")
    d << gf.components.circle(layer=LAYER.B_Cu, radius=5)
    dref = c << d
    dref.move(origin=(0, 0), destination=(4, 0))
    dref2 = c << d

    # c = gf.components.text(layer=LAYER.F_Cu)
    # c = LAYER_VIEWS.preview_layerset()

    gerber_path = PATH.repo / "extra" / "gerber"
    gerber_path.mkdir(exist_ok=True, parents=True)

    # This requires that the PCB technology (commented-out code above) is installed
    c.show()

    to_gerber(
        c,
        dirpath=gerber_path,
        layermap_to_gerber_layer=layermap_to_gerber,
        options=GerberOptions(resolution=1e-6),
    )
